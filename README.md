> [!NOTE]
> (From CoCoA notes)
> `v4.11.1` benchmark: do not include CAMB (or the Hybrid Emulator); **Includes TATT in** ($\xi_{\pm}, \gamma_t$) **and non-limber in** $w_{gg}(\theta)$.
> CLOE-LIB caveat: We were not able to make TATT work on CloeLib. We were also not able to speed-up cloelib with `OMP_NUM_THREADS` flag
>
> CPU: `Intel(R) Core(TM) i9-10940X CPU @ 3.30GHz` (`1/8 OpenMP cores`). *Times are approximate*.
>
> - **LSST-Y1-Real 3x2pt**: (CoCoA) `0.29/0.06s`, (DESC-CCL)`7.96/1.72s`, (CLOE-LIB) 0.23/0.23s. **CoCoA speed-up (CCL)**: `27/28x`
> - **Roman-Real 3x2pt**: (CoCoA) `0.45/0.095s`, (DESC-CCL) `8.17/1.96s`, (CLOE-LIB) 0.27/0.27s. **CoCoA speed-up (CCL)**: `18/20x`
> - **Roman-Fourier 3x2pt**:  (CoCoA) `0.08/0.03s`, (DESC-CCL) `0.65/0.36s`. **CoCoA speed-up**: `7.5/21x`
> - **DES-Y3xPlanck 6x2pt**  (CoCoA) `0.40/0.075s`
> - **DES-Y3-Real 3x2pt (des_y3 repo)**  (CoCoA)~`0.25/0.05s`
> 
