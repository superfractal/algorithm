# v199 local Release build

Created by GPT-6 on 2026-09-12.

Applied the completed v199 GPU locator implementation and built the production GUI with Release optimization, as explicitly requested by the user. All four GPU optimizations are enabled by default. Existing application, batch enumeration, rendering defaults, title and timing UI were preserved.

Build completed using Clang and Ninja, with `-O3 -DNDEBUG`; locator and GPU arithmetic retain `-fno-fast-math -ffp-contract=off -fopenmp`. Eight compiled orbit shader variants are included in `bin/` and the locator source directory. Configuration used the existing local nativefiledialog source without fetching dependencies. OpenCV headers emitted NaN/infinity warnings under the existing renderer fast-math flags; the build succeeded.

Executable: `Release/bin/RFF.exe` (4493824 bytes).

SHA-256: `e385cb4e6e8eff40c64f1eb7d82eb1b0a7962c7a3defb845ff94031ee4716f8a`.

Complete pre-application Release backup: `Debug/apply_v199_001/Release_before/`; its files were hash-checked before application. Build log and metadata: `Debug/apply_v199_001/build.log` and `BUILD.json`.

No executable was launched and no tests or performance measurements were run for this Release build, per the user's instruction. Previous Debug experiment evidence is in `Debug/source/v199/RESULT.md`; it is not validation of this newly built GUI executable. This build does not establish an all-input speedup or achievement of the research targets. No commit or push was performed.
