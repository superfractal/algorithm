<!-- Created by GPT-6 on 2026-09-12 -->
<!-- Modified by GPT-6 on 2026-09-12 -->
# v202 Release application and evidence

Applied the tested, Release-optimized v202 executable, its matching project helper DLL, the v202 locator sources, and reviewed source/reference documentation to Release. The user waived remaining checks on 2026-09-12 after stopping Computer Use. No additional GUI checks or benchmarks were run after that instruction. No commit or push was performed.

## Completed measurements

The before executable was the local v201 Release, not the older dd8f1c3 commit. The candidate was built in the complete production-equivalent copy `Debug/source/real/v202` with Clang/Ninja, Release `-O3 -DNDEBUG`; strict locator arithmetic retains `-fno-fast-math -ffp-contract=off`. The exact tested executable and project helper DLL were installed, without rebuilding different binaries afterward.

Hardware: AMD Ryzen 7 9700X, NVIDIA GeForce RTX 5060 Ti, Vulkan driver 2559967232. Both GUI locator runs used 15 threads, 1280x720, the CPU image renderer, and the same default settings. Locator GPU dispatch remains separately automatic for eligible deep inputs. GUI logs do not record which individual stages used GPU or fallback.

Per-file locator seconds exclude FPG/reference preparation, rendering, and saving. Baseline: `Debug/locate_minibrot_test/batch_output_0016`; candidate: `batch_output_0017`. Both completed all 16 supported inputs and ended `done.`. Automatic FPG was computed from original locations without Move To Center or reference reuse.

| Input | Before seconds | v202 seconds | Speedup |
| --- | ---: | ---: | ---: |
| 1.rfl | 0.259654100 | 0.192189600 | 1.351x |
| 2.rfl | 4.422613200 | 3.612697800 | 1.224x |
| 3.rfl | 0.046525100 | 0.037261100 | 1.249x |
| 4.rfl | 6.023915900 | 3.873582900 | 1.555x |
| 5.rfl | 0.001523700 | 0.000867600 | 1.756x |
| 6.rfl | 0.000891400 | 0.000896900 | 0.994x |
| 7.rfl | 2.951114200 | 2.137611100 | 1.381x |
| 8.rfl | 1.462279100 | 1.137435300 | 1.286x |
| 9.rfl | 0.746647500 | 0.576612200 | 1.295x |
| 10.rfl | 4.308548900 | 3.644439900 | 1.182x |
| 11.rfl | 12.879560200 | 9.274184100 | 1.389x |
| 12.rfl | 29.165584900 | 27.717353100 | 1.052x |
| 13.rfl | 112.078507000 | 112.502581300 | 0.996x |
| 14.rfl | 0.478870100 | 0.441400400 | 1.085x |
| 16.rfl | 27.138509600 | 17.331384400 | 1.566x |
| 17.rfl | 7.275807900 | 6.538319800 | 1.113x |

Inputs 1-4 improved in this single GUI pair. Inputs 5/6 have no speed gate. Input 13 was 0.38% slower in this sample; this is not evidence of universal acceleration. Input 17 met its ten-second target. The higher 10x/2x research goals are not achieved by these results.

All 16 saved RFL files are byte-for-byte identical to the baseline. Their log coordinates and float32 zooms match their saved files. Prior v202 bounded numerical evidence is retained in `Debug/source/v202/RESULT.md`; it is not relabeled as fresh testing. There is no independently supplied expected location for the later fixtures.

Fresh visual inspection showed minibrots for saved outputs 1-5. Screenshots retained in `Debug/apply_v202_001/visual/`. Case 6 and additional fresh deep-output visual checks were not completed after the user stopped Computer Use; the user subsequently waived remaining checks. This is an application based on completed numerical/render regression evidence and the user's reduced verification scope, not a claim that every original four-stage check was completed.

SHA rendering regression: all 17 inputs completed in fresh baseline `Debug/sha256_test/batch_output_0013` and candidate `batch_output_0014`. All RFM/PNG hashes match their own logs and each other: zero differing input entries. Settings match. The historical `batch_output_0001` was cancelled and was not used as a valid baseline.

Input 15 remains outside the known period limit; input 18 was outside this v202 test scope. Both were temporarily held during locator batches. All 18 original location files were restored with their original SHA-256 values, recorded in `fixture_restore.log`.

## Source and reference review

Review scope: local changes since `dd8f1c32bad8cc4429b0c531b27d58731755c849`, including untracked GPU implementation files. [SOURCES_AND_REFERENCES.md](SOURCES_AND_REFERENCES.md) now records the attachment, local v197-v202 implementation lineage, dependency APIs, and generated-shader inputs. [About_the_Algorithm.md](About_the_Algorithm.md) describes current v202 behavior. Missing historical-document pointers and misleading inherited prototype comments were corrected. Existing author history, license bodies, and the historical attachment document were preserved.

The reviewed implementation lineage and inherited references are recorded in [SOURCES_AND_REFERENCES.md](SOURCES_AND_REFERENCES.md) and the [consolidated bibliography](About_the_Algorithm.md#10-consolidated-references). The earlier detailed audit is retained in the local archive.

## Installed artifact identity and recovery

- `bin/RFF.exe` SHA-256: `0786440d0e3bfb7421afd999acbd0bb466b1cb4b2d90619ba66bdaa4ea07d47f`
- `bin/libvulkan_helper.dll` SHA-256: `b3bf3a7f2fb14c047e1fc49e9b0f9dfc3a72261d953cb7ffe3a6a416e0bf09ff`
- Complete before backup: `Debug/apply_v202_001/Release_before/` (2512 files), with `backup_manifest.json`.
- Exact promoted-file hashes: `Debug/apply_v202_001/PROMOTED.json`.
- Source identity and build inputs: `STAGED_BUILD.json`, `PREPARE.json`, `build_stage.log`, `build_stage_documented.log` in that evidence directory.
- Detailed comparisons: `locator_comparison.json`, `sha_comparison.json`, `NUMERICAL_REVIEW.md`, `SOURCE_AUDIT.md`.

The helper DLL was rebuilt from unchanged project source as part of the tested GUI build and installed with its matching executable. Other DLL/SPIR-V dependencies match the existing Release. Unrelated user changes, including documentation deletions, were preserved. Installed hashes were checked; the installed-path application was not relaunched after the user's instruction to omit remaining checks.
