# v201: input 17 performance and precision guard

Created by GPT-6 on 2026-09-12.

The original shallow/long-period dispatch forced return jets even where they were expensive. v201 limits individual jet build/replay attempts to a polled 1.5-second dispatch budget, restores the original parameter and discards partial maps on expiration, then uses the direct multiple-shooting proposal. Full bounded verification remains mandatory. Input 4 retains cached return jets. A separate guard now clamps precision estimates before integer conversion, preventing undefined conversion on a degenerate zero-orbit trial without reducing ordinary precision.

Input 17: 7.2404641 seconds from the original input with automatic FPG, full bounded verification and an actual CPU iteration map showing a minibrot. Three additional final-version runs: 7.1191316, 7.0868552, 7.0763992 seconds. All meet the requested <=10 seconds. These are optimized headless locator-only timings, excluding FPG/reference/MPA preparation, rendering and saving; the rebuilt GUI was not timed. The old run's initial jet build alone measured 79.432 seconds, but the old locator did not complete before it was stopped, so no precise complete-run speedup ratio is claimed.

Fresh comparisons on inputs 1-14 and 16 exactly match v200 candidate coordinates, period, residual and bounded zoom endpoints. Additional checks: 120 automatically generated known-nucleus probes and 163860 precision-estimate parity cases with sanitizer edge checks. Existing-input renders and the GUI batch/SHA-256 suite were not repeated; input 17 was fully rendered in the optimized harness. See `Debug/source/v201/RESULT.md`, `VERIFICATION.json` and `INPUT17_MASK.png` for details and scope.

Input 15 remains unsupported because FPG period 126294672 exceeds the unchanged 100000000 limit. This is not an all-location correctness guarantee or an all-input speedup claim. Generated solver-test failures are saved as RFL files by `Debug/robustness_v200_001/probe.cpp` for reproducibility.

Release executable: `Release/bin/RFF.exe` (4516864 bytes), SHA-256 `9706363b9261f00225bc0a95d1711cc258f03a65aa05c8562877d17921a1a71e`. Release optimization and strict locator arithmetic flags were confirmed. The CMake cache directory remains named `build_v199`; its objects now contain v201. This task did not launch the rebuilt GUI. Application/rendering code, defaults, title and normal batch enumeration were preserved.

Full pre-application backup: `Debug/apply_v201_001/Release_before/`. Build metadata and log: `Debug/apply_v201_001/BUILD.json` and `build.log`. No web search, sub-agents, commit or push were used.
