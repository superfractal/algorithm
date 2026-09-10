# STMS-RJ: sensitivity-tapered multiple shooting + cached return jets

This is an experimental **nucleus-finding core**, not a complete replacement for RFF's locate-minibrot UI operation. In particular, the RFF-compatible output zoom is not implemented. The empirical tests supplied with this package do not constitute a proof of global convergence or a same-machine comparison against a release binary.

## 1. Problem and period candidate

For `f_c(z)=z*z+c`, evolve the critical orbit `z_0=0` and its parameter derivative:

```
z_(n+1) = z_n² + c
D_(n+1) = 2 z_n D_n + 1,  D_0 = 0
```

The first index satisfying `|z_n/D_n| < 32 * 10^(-input_log_zoom)` is a **candidate** period. This search radius is a heuristic, not a guarantee that an arbitrary viewport has a unique desired minibrot. The supplied six periods are discovered from the inputs; they are not present as constants or looked up by filename.

After selecting `p`, Newton's correction is `dc = -z_p/D_p`. The root sought is the critical-period nucleus, not necessarily the partially refined center returned by another program's stopping convention.

## 2. Sensitivity-weighted fixed-point precision

Let `A_n = product(2*z_j, j=1..n-1)`. Early orbit errors are amplified roughly in proportion to `A_n`. The fixed-point implementation chooses the number of fractional bits from

```
fractional_bits(n) ≈ goal_bits - log2|A_n| + guard_bits
```

rounded to a 64-bit limb boundary and clipped between a floor of 128 bits and the allocated maximum. Large integer exponents are separated from a double mantissa for cheap pilot derivatives. Orbit values themselves are **not** stored as ordinary doubles.

This replaces the cost of uniformly evaluating every orbit step at the final center precision. If the coordinate is genuinely very close to -2, the real state is additionally represented as `integer_part + small_residual`, and the square is expanded algebraically. That optimization is selected from the coordinate, not from an input number.

The guard margins and double-valued logarithmic estimates are empirical engineering choices. They are not outward-rounded interval arithmetic. A production integration should retain a verified fallback when the reference orbit becomes ill-conditioned.

## 3. Parallel multiple shooting

The pilot saves a state and sensitivity every 256 steps. Divide the orbit into blocks. At each Newton step each block independently evaluates

```
f_i(s_i,c), A_i = ∂f_i/∂s_i, B_i = ∂f_i/∂c
r_i = f_i(s_i,c) - s_(i+1)
```

with fixed endpoint states `s_0=s_end=0`. Parallel block evaluations are followed by the serial prefix recurrences

```
t_0 = 0, u_0 = 0
t_(i+1) = A_i t_i + r_i
u_(i+1) = A_i u_i + B_i
dc = -t_end/u_end
ds_i = t_i + u_i dc
```

Both the parameter and intermediate states are corrected. That distinction is essential: parallel propagation of independent chunks without correcting their boundary states would not solve the same orbit equation.

The full-precision coordinate is allocated once. Orbit precision and derivative precision increase as the Newton residual shrinks. A final evaluation with inexpensive exponent-extended derivatives is used when the predicted remaining correction is already small. The stopping condition checks the parameter correction **and** sensitivity-normalized state corrections.

Default parallelism is four OpenMP threads. Default block counts are eight per thread below 100,000 steps and twelve per thread above it, limited by available checkpoints.

## 4. Long low-precision orbits: Taylor return jets

When input_log_zoom is below 500 and no period candidate appears in the first 100,000 pilot steps, use cached return jets. The dispatch pilot remains inside the measured time.

The first state dependence is naturally even. Write `w=z_initial²` and represent a return map as

```
S_p(w) = f_c^p(z_initial)
B_p(w) = ∂S_p(w)/∂c
```

truncated at degree 24 in `w`. Initially `S_1=c+w` and `B_1=1`. A micro-step updates

```
S <- S²+c
B <- 2*S_old*B+1
```

When the orbit makes a new record-small return, save its map. To jump by a cached return map `M`:

```
W = S²
S_new = M.S(W)
B_new = 2*S*M.S'(W)*B + M.B(W)
```

This is a nonlinear return map, not just a linear perturbation table.

### Domain and truncation estimates

The implementation tracks a positive rational majorant of the form

```
deviation(r) <= G*r / (1-H*r)
```

and a domain cap preventing internal reference returns from reaching zero. Before a micro-step, restrict the domain using `eta=1/64`; then update

```
G_new = 2*|z|*G + 1
H_new = H + G²/G_new
```

A saved map records a domain radius `R` and amplitude majorant `L`. For `q=|w|/R`, its state Taylor tail is estimated by

```
L*q^(K+1)/(1-q), K=24.
```

Jump acceptance normalizes that error by the global parameter sensitivity. During period discovery the test is additionally restricted by a parameter disk of radius `32*10^-input_log_zoom`, to reduce the risk of skipping a candidate return. These tests are **not interval-certified**: polynomial roundoff, accumulated jump error, and floating-point log estimates are not enclosed by a formal proof. The delivered validation measures the actual accuracy on the supplied examples.

## 5. Reuse maps across Newton steps

The substantial case-4 improvement comes from **not rebuilding the jets at every Newton step**. Save maps at a reference parameter `c0`. For the new parameter `c=c0+dc`, approximate

```
M_c(w) = M_c0(w) + dc * ∂M_c0(w)/∂c.
```

Use `s=|dc|/R` and reject a map unless the estimated omitted state and parameter terms satisfy

```
error <= L * (q^(K+1)+s²) / ((1-q)*(1-s)).
```

An inaccurate large map is not forced through: the evaluator tries a smaller map, and otherwise takes a direct micro-step. In the case-4 development run, later Newton evaluations needed only 37–365 scalar map jumps and no direct micro-steps.

This reuse is the difference between a roughly 20-second prototype that rebuilt a degree-20 table on every correction and the final roughly 5-second center-search prototype. Those development figures are not same-source controlled ablation benchmarks; the repeatable final measurements are in `results/`.

## 6. Size information, not a display-zoom replacement

At the nucleus, the result includes

```
A = product(2*z_j, j=1..p-1)
D = dz_p/dc
log10|A*D| and arg(A*D).
```

The local quadratic return-map scaling motivates using `C_effective = A*D*delta_c`. This is useful intrinsic size information, but the display logZoom also depends on the program's viewport conventions, boundary/escape test, finite iteration cap, and stopping rule. `output_log_zoom` is therefore deliberately `null`; an arbitrary fitted constant from the expected files is not substituted for those missing conventions.

## 7. Novelty and reference boundary

Newton iteration, multiple shooting, sensitivity propagation, Taylor composition, and Cauchy/majorant estimates are established mathematical ideas. This package claims an experimentally useful **combination and implementation**, not the invention or priority of those general methods. Its formulas were derived for this task. No unconsulted paper is presented as a source that was actually used.

The project sources consulted for the existing implementation and build context were:

- Merutilm/RFF-2.0 repository README: https://github.com/Merutilm/RFF-2.0
- `src/rff2/mb/MB2Locator.cpp` on the master branch: https://raw.githubusercontent.com/Merutilm/RFF-2.0/master/src/rff2/mb/MB2Locator.cpp
- `CMakeLists.txt` on the master branch: https://raw.githubusercontent.com/Merutilm/RFF-2.0/master/CMakeLists.txt

The repository checkout and the dependencies needed for the GUI release could not be acquired in the execution environment. The consulted master source was not pinned to the exact release used for the supplied timings. The new solver is independently written and has not been compiled inside RFF.
