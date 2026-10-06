#!/bin/bash
# Leave this private environment, keeping the Conda base active.
if [[ "${BASH_SOURCE[0]}" == "$0" ]]; then
  echo "Use: source stop_ccl.sh" >&2
  exit 1
fi
if [[ "${VIRTUAL_ENV:-}" != \
      "$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)/.local" ]]; then
  echo "This benchmark's private environment is not active." >&2
  return 1
fi
deactivate || return 1
return 0
