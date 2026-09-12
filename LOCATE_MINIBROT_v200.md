# v200: fixes for inputs 14 and 16

Created by GPT-6 on 2026-09-12.

Fixed a zero-length pilot partition when a period is divisible by 256, and added bounded trial-radius refinement with full re-verification for slightly insufficient connection budgets. Full precision, connection containment, residual and zoom validation remain mandatory. Only three locator files changed; rendering, UI, defaults and batch enumeration were preserved.

Optimized headless validation from the original RFL files: 14 completed in 0.425246 seconds, 16 in 26.0551 seconds, with visible minibrots in actual CPU iteration maps. Inputs 1-13 produced byte-identical saved locations and 1280x720 iteration maps against v199. These are locator-only source timings, not measurements of this rebuilt GUI executable. Numerical boundary, rejection, directed enclosure, CPU/GPU and cancellation checks are documented in `Debug/source/v200/RESULT.md` and `VERIFICATION.json`.

Additional input 15 remains unresolved: automatic FPG period 126294672 exceeds the existing 100000000 limit. The 16-entry batch completed but did not pass all 16. No independently supplied expected result exists for 14 or 16. No new all-input speedup claim is made.

Production executable: `Release/bin/RFF.exe`, 4511744 bytes, SHA-256 `17adf6ae5dc25daf00576afa38f75b1f69109f07ec79322fc8c382d0de10b7f7`. Built with Release optimization; existing `build_v199/` is the CMake build-directory name and now contains the updated v200 objects. The GUI was not launched, and GUI batch/SHA-256 suite checks were not run in this repair. The full iteration-map regression above was run in the optimized harness.

Complete pre-application Release backup: `Debug/apply_v200_001/Release_before/`. Build log and metadata: `Debug/apply_v200_001/build.log` and `BUILD.json`. No commit or push was performed.
