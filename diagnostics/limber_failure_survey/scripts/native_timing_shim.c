/* Isolated macOS diagnostic: retain CCL's native C integrand while switching
 * only an explicitly armed, already-warmed Limber call from QAG to CQUAD.
 * Never preload this library into a production process. No Python callbacks.
 * Mode 0: transparent. Mode 1: native QAG + native fallback, recorded.
 * Mode 2: CQUAD with a private reusable workspace per OpenMP thread.
 * CCL's outer QAG workspace is still allocated in mode 2; public-call timing
 * includes this unused allocation. Kernel timing excludes workspace allocation.
 */
#include <dlfcn.h>
#include <stdatomic.h>
#include <stdint.h>
#include <stddef.h>
#include <time.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_integration.h>

typedef int (*qag_fn)(const gsl_function *, double, double, double, double,
                     size_t, int, gsl_integration_workspace *, double *, double *);
typedef int (*cquad_fn)(const gsl_function *, double, double, double, double,
                       gsl_integration_cquad_workspace *, double *, double *, size_t *);
static qag_fn original_qag;
static cquad_fn original_cquad;
static _Atomic int armed_mode;
static _Atomic uintptr_t expected_callback;
/* Written only before/after a synchronous angular_cl call; never during it. */
static double relative_tolerance_override;
static _Thread_local gsl_integration_cquad_workspace *thread_cquad;
static _Thread_local size_t thread_capacity;
enum {QAG_CALLS, CQUAD_CALLS, UNEXPECTED_CALLBACKS, WORKSPACE_ALLOCATIONS,
      KERNEL_NANOSECONDS, QAG_FAILURES, CQUAD_FAILURES, FALLBACK_CQUAD_CALLS,
      COUNTER_COUNT};
static _Atomic uint64_t counters[COUNTER_COUNT];
static _Atomic uint64_t status_counts[2][64];

static uint64_t nanoseconds(void) {
  struct timespec t;
  clock_gettime(CLOCK_MONOTONIC, &t);
  return (uint64_t)t.tv_sec * UINT64_C(1000000000) + (uint64_t)t.tv_nsec;
}

static int selected_callback(const gsl_function *f) {
  uintptr_t observed = (uintptr_t)f->function;
  uintptr_t expected = atomic_load(&expected_callback);
  if (!expected) {
    atomic_compare_exchange_strong(&expected_callback, &expected, observed);
    expected = atomic_load(&expected_callback);
  }
  if (observed == expected) return 1;
  atomic_fetch_add(&counters[UNEXPECTED_CALLBACKS], 1);
  return 0;
}

static void record_call(int method, int status, uint64_t elapsed) {
  atomic_fetch_add(&counters[method == 0 ? QAG_CALLS : CQUAD_CALLS], 1);
  atomic_fetch_add(&counters[KERNEL_NANOSECONDS], elapsed);
  if (status != GSL_SUCCESS)
    atomic_fetch_add(&counters[method == 0 ? QAG_FAILURES : CQUAD_FAILURES], 1);
  if (status >= 0 && status < 64)
    atomic_fetch_add(&status_counts[method][status], 1);
}

static int replacement_qag(const gsl_function *f, double a, double b,
                           double epsabs, double epsrel, size_t limit, int key,
                           gsl_integration_workspace *w,
                           double *result, double *error) {
  int mode = atomic_load(&armed_mode);
  if (!mode || !selected_callback(f))
    return original_qag(f, a, b, epsabs, epsrel, limit, key, w, result, error);
  if (mode == 1) {
    uint64_t start = nanoseconds();
    int status = original_qag(f, a, b, epsabs, epsrel, limit, key, w, result, error);
    record_call(0, status, nanoseconds() - start);
    return status;
  }
  if (!thread_cquad || thread_capacity < limit) {
    if (thread_cquad) gsl_integration_cquad_workspace_free(thread_cquad);
    thread_cquad = gsl_integration_cquad_workspace_alloc(limit);
    thread_capacity = thread_cquad ? limit : 0;
    atomic_fetch_add(&counters[WORKSPACE_ALLOCATIONS], 1);
  }
  if (!thread_cquad) {
    record_call(1, GSL_ENOMEM, 0);
    return GSL_ENOMEM;
  }
  size_t evaluations = 0;
  double tolerance = relative_tolerance_override > 0 ? relative_tolerance_override : epsrel;
  uint64_t start = nanoseconds();
  int status = original_cquad(f, a, b, epsabs, tolerance, thread_cquad,
                              result, error, &evaluations);
  record_call(1, status, nanoseconds() - start);
  return status;
}

static int replacement_cquad(const gsl_function *f, double a, double b,
                             double epsabs, double epsrel,
                             gsl_integration_cquad_workspace *w,
                             double *result, double *error, size_t *evaluations) {
  int mode = atomic_load(&armed_mode);
  if (!mode || !selected_callback(f))
    return original_cquad(f, a, b, epsabs, epsrel, w, result, error, evaluations);
  /* This path is CCL's own QAG-roundoff fallback. Direct mode 2 above calls
   * original_cquad, so it cannot be counted twice. */
  uint64_t start = nanoseconds();
  int status = original_cquad(f, a, b, epsabs, epsrel, w, result, error, evaluations);
  atomic_fetch_add(&counters[FALLBACK_CQUAD_CALLS], 1);
  record_call(1, status, nanoseconds() - start);
  return status;
}

#define INTERPOSE(replacement, replacee) \
  __attribute__((used)) static struct { const void *replacement; const void *replacee; } \
  interpose_##replacee __attribute__((section("__DATA,__interpose"))) = \
  { (const void *)(uintptr_t)&replacement, (const void *)(uintptr_t)&replacee }
INTERPOSE(replacement_qag, gsl_integration_qag);
INTERPOSE(replacement_cquad, gsl_integration_cquad);

__attribute__((constructor)) static void resolve_originals(void) {
  /* dyld excludes the interposing image's own references from replacement.
   * RTLD_NEXT lookup can itself resolve the interposed symbol on macOS. */
  original_qag = &gsl_integration_qag;
  original_cquad = &gsl_integration_cquad;
}

int limber_timing_ready(void) {
  return original_qag && original_cquad && original_qag != replacement_qag
         && original_cquad != replacement_cquad;
}

/* Caller must have joined all previous angular_cl work before changing mode.
 * Guards and cached workspaces persist across repetitions; counters do not. */
int limber_timing_begin(int mode, double epsrel) {
  if (!limber_timing_ready() || mode < 1 || mode > 2 || atomic_load(&armed_mode)) return -1;
  for (int i = 0; i < COUNTER_COUNT; ++i) atomic_store(&counters[i], 0);
  for (int i = 0; i < 2; ++i)
    for (int j = 0; j < 64; ++j) atomic_store(&status_counts[i][j], 0);
  relative_tolerance_override = epsrel;
  atomic_store(&armed_mode, mode);
  return 0;
}
void limber_timing_end(void) { atomic_store(&armed_mode, 0); }
uint64_t limber_timing_counter(int index) {
  return index >= 0 && index < COUNTER_COUNT ? atomic_load(&counters[index]) : 0;
}
uint64_t limber_timing_status_count(int method, int status) {
  return method >= 0 && method < 2 && status >= 0 && status < 64
    ? atomic_load(&status_counts[method][status]) : 0;
}
/* Only frees the calling thread's cache. Worker caches are bounded by the
 * configured thread count and intentionally live until this diagnostic exits. */
void limber_timing_release_current_thread(void) {
  if (thread_cquad) gsl_integration_cquad_workspace_free(thread_cquad);
  thread_cquad = NULL;
  thread_capacity = 0;
}
