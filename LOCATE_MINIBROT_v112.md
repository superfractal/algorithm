# Locate Minibrot v112: STMS-RJ integration

Created by GPT-6 on 2026-09-10.

PASS and locally applied under the user's revised criterion: the intended
minibrot is visible on all six inputs, and 1/2/4 are faster than current Release.
Input 3 is correctness-only for this request; 5/6 are shallow stability tests.
The older 10x/2x research target is separate and is not claimed achieved.

## Measured GUI performance

Original Release: Debug/locate_minibrot_test/batch_output_0003.
Final optimized v112: batch_output_0005. Both complete with done.
The timer covers locateMinibrot, including full-precision validation and its
internal result reference/table construction. Input preparation, rendering
and saving are excluded. These are same-machine single-run measurements,
not results from an unoptimized Debug build or from another fractal program.

| Input | Original seconds | v112 seconds | Speedup | Role |
| --- | --- | --- | --- | --- |
| 1 | 11.367341200 | 1.708921400 | 6.652x | Speed + correctness |
| 2 | 130.797897600 | 20.082633900 | 6.513x | Speed + correctness |
| 3 | 0.163850100 | 0.116658200 | 1.405x | Correctness/stability |
| 4 | 57.437732400 | 21.264319100 | 2.701x | Speed + correctness |
| 5 | 0.002756300 | 0.003419600 | 0.806x | Correctness/stability |
| 6 | 0.002955200 | 0.003267500 | 0.904x | Correctness/stability |

The first v112 GUI run (batch_output_0004) also passed. After adding graceful
handling of rejected proposals, the final build produced byte-identical RFLs
on all six inputs. Those are the v112-to-v112 hashes, not equality to the old
locator's coordinates. Original and new centers need not match digit-for-digit.

## Correctness and regression

All six v112 outputs were opened through the GUI and rendered to completion;
the intended centered minibrot was visible in each. Twelve independent GMP
diagnostics (headless and GUI saved files) passed full-orbit, derivative,
critical-product and all-proper-divisor checks. Saved float zooms equal the
independent computed zoom rounded to float. The runtime MPFR acceptance pass
uses uniform full precision without jets/tapering and requires |A Fp|<1e-12.
These numerical checks are not interval root proofs.

Rendering SHA-256 regression: all 17 fixed inputs completed; every RFM/PNG
artifact matches its own logged hash; zero differing inputs between original
Release batch_output_0002 and v112 batch_output_0003 under Debug/sha256_test.
The pre-existing output 0001 was cancelled after three inputs and is preserved;
it could not be used as a complete baseline. Fresh output 0002 replaces it for
this comparison. The comparison uses identical fixed SHA input positions,
not the different positions produced by the two locator algorithms.

MPA source files are identical between original Release and this candidate.
Only the locator and its numerical support were changed. The tested executable
and rebuilt Vulkan helper DLL are deployed together. Rendering regressions
include both fresh-reference and reference-reuse cases.

Environment: NVIDIA GeForce RTX 5060 Ti, driver 2559967232, 15 threads, CPU mode,
1280x720. Locator: fresh defaults, reference reuse disabled. SHA: Classic 1,
noise reduction off, linear palette interpolation off, reference/MPA
compression off, alternating reference reuse as defined by the batch runner.
Compiler: Clang 22.1.4; CMAKE_BUILD_TYPE=Release. Numerical bridge flags end in
-fno-fast-math -ffp-contract=off -fopenmp while retaining -O3 -DNDEBUG
-march=native. Experiment flags and source provenance are preserved under
Debug/source/v112. No subagents or web search were used.

## Implementation and limits

The attached sensitivity-tapered multiple shooting and cached degree-24 Taylor
return jets propose the nucleus. The discovered period must match RFF FPG.
The independent MPFR pass supplies logZoom=log10|A D|+2; the adapter preserves
the GUI's existing 1.5 framing-offset convention. No expected center or
input-number-specific branch is used. Failed proposals return nullptr and log
the reason. A mutex isolates the persistent map cache/precision setting;
long loops poll cancellation and a ten-minute deadline. Safety tests cover
exact nuclei, period mismatch rejection, cancellation and post-cancel cleanup.

The prototype's error estimates are heuristic; the full-precision acceptance
pass is mandatory. The supplied prototype supports input logZoom 0..20000
and period at most 100 million; out-of-range inputs or failure to converge
report failure. The existing input bounds and fallback behavior should be
expanded only with additional correctness tests, not by accepting weaker
residuals. No legacy Newton fallback was silently substituted in these tests.

Build with Release configuration. MPFR prefers the static libmpfr.a; the
tested OpenMP runtime libomp.dll is included in bin. The numerical bridge
requires OpenMP-capable Clang and GMP/MPFR. The complete pre-change Release
backup is Debug/apply_v112_001/Release_before, with BACKUP_MANIFEST.txt.
No commit or push was performed under the local application authorization.

## Sources actually used

- User attachment locate_minibrot_experiment.zip, SHA-256
  25a4d9bde0496ffdb6f5e3ba7c79cfaa703a8035a885fe7fd001e56414585750.
  Its original algorithm description is copied as
  [LOCATE_MINIBROT_STMS_RJ_ALGORITHM.md](LOCATE_MINIBROT_STMS_RJ_ALGORITHM.md).
  The full unchanged archive contents are under Debug/source/v112/attachment.
- Local RFF MB2Locator.cpp, MB2RenderData.hpp, FnExplore.cpp and RFF2.cpp
  for FPG, precision, framing and GUI integration.
- Local v81 quality.hpp/inverse_samples.hpp/interval_solver.hpp for the
  independent uniform-precision acceptance pass and size calculation.
- Local v76 audit_gmp.cpp for separate saved-root diagnostics.

This is an adaptation of the supplied algorithm, not a claim of novelty or
priority for Newton iteration, multiple shooting, Taylor maps or majorants.
The upstream links in the supplied algorithm document are its author's source
history; they were not browsed during this implementation.

## Deployed hashes

- bin/RFF.exe: `19cc0d3f5f4ac022f13d5bf3142eecd9c3681c89468e4ab242c6dab0bc766215`
- bin/libvulkan_helper.dll: `4c847e3783c7c99a7505f2f50bf02b7918c05a548b6f9ea9e4b1ce6471cd5cd6`
- bin/libomp.dll: `f9ac1881ed57f9cf6277b55622059afcce55626495821651d83f2cd8d62f4dd7`
