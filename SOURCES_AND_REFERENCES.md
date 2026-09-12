<!-- Created by GPT-6 on 2026-09-12 -->
<!-- Modified by GPT-6 on 2026-09-12 -->
# Locate Minibrot source index

The complete bibliography is in [About the Algorithm, section 10](About_the_Algorithm.md#10-consolidated-references). That document uses the [superfractal/algorithm reference at a pinned commit](https://github.com/superfractal/algorithm/tree/90a006d45e88b204b0ab9b8f5d0c18d7ea928a79) and describes the current v202 additions.

## Original proposal

The user-supplied `locate_minibrot_experiment.zip` was adapted in v112 on 2026-09-10. Its SHA-256 is `25a4d9bde0496ffdb6f5e3ba7c79cfaa703a8035a885fe7fd001e56414585750`. The adapted files are `prototype/locate.cpp`, `locate.hpp`, `numerics.inc` and `return_jets.inc` under `src/rff2/mb/`.

The [original description](LOCATE_MINIBROT_STMS_RJ_ALGORITHM.md) and [v112 integration record](LOCATE_MINIBROT_v112.md) retain that history. The original ZIP is not distributed. Current differences include degree-20 jets, automatic FPG period, bounded verification and output zoom; see About the Algorithm for the active call path.

## Local additions

Historical `Debug/source/vN/` records named below are now stored under `_local_only/Debug/source/vN/` in the development checkout. The Release version reports remain alongside this document.

| Implementation | Local source and scope |
| --- | --- |
| [integer_square.hpp](src/rff2/mb/prototype/integer_square.hpp) | The user's exact two-product complex-square proposal, implemented and calibrated in `Debug/source/v197/DESIGN.md`. It conditionally uses `(x+y)(x-y)` and `xy` with unchanged caller truncation points. This is established algebra, not a newly claimed identity. |
| [ntt.comp](src/rff2/mb/ntt.comp), [ntt_gpu.hpp](src/rff2/mb/ntt_gpu.hpp), [ntt_gpu.cpp](src/rff2/mb/ntt_gpu.cpp) | Local v198 integer NTT implementation, using v197's Vulkan device/buffer setup and the user's GPU direction. The local record is `Debug/source/v198/DESIGN.md`. Two-prime NTT, Montgomery arithmetic, exact CRT reconstruction, and carry scans are established techniques. No newly consulted external paper or copied NTT library was identified in the inspected development record. |
| [orbit_gpu.hpp](src/rff2/mb/orbit_gpu.hpp), [orbit_gpu.cpp](src/rff2/mb/orbit_gpu.cpp), [orbit_phases.inc](src/rff2/mb/orbit_phases.inc), [gpu_search.inc](src/rff2/mb/prototype/gpu_search.inc), [gpu_segments.inc](src/rff2/mb/verified_quality/gpu_segments.inc) | Local v198 fixed-point orbit integration and bounded CPU summaries, followed by v199 optimization. The orbit retains full fixed-point state; 128-bit summaries carry an explicit omitted-tail bound for derivative processing. Sources/evidence: v198 DESIGN/RESULT and `Debug/source/v199/RESULT.md`; the [v199 Release record](LOCATE_MINIBROT_v199.md) has its own build-only scope. |
| [optimized_phases.inc](src/rff2/mb/optimized_phases.inc) and [make_orbit_shader.py](src/rff2/mb/make_orbit_shader.py) | Local v199 work implementing the user's four requested GPU experiments: reductions, combined transform-space outputs, fused radix-4 stages, and overlapping GPU submission with CPU bound updates. `Debug/source/v199/RESULT.md` identifies the v198 implementation as its source and does not claim new external research sources. |
| [STMSBridge.cpp](src/rff2/mb/STMSBridge.cpp), [numerics.inc](src/rff2/mb/prototype/numerics.inc), [segments.hpp](src/rff2/mb/verified_quality/segments.hpp) — v200 changes | Local diagnosis of terminal checkpoint ordering and insufficient trial connection radii. Refinements trigger recomputation without relaxing acceptance. Sources/evidence: `Debug/source/v200/RESULT.md` and [LOCATE_MINIBROT_v200.md](LOCATE_MINIBROT_v200.md). |
| [locate.cpp](src/rff2/mb/prototype/locate.cpp), [return_jets.inc](src/rff2/mb/prototype/return_jets.inc), [numerics.inc](src/rff2/mb/prototype/numerics.inc) — v201 changes | Local workload profiling and precision-conversion checks. A bounded jet attempt can restore the original state and fall back to the direct proposal; precision estimates are bounded before narrowing. Sources/evidence: `Debug/source/v201/RESULT.md` and [LOCATE_MINIBROT_v201.md](LOCATE_MINIBROT_v201.md). |
| [compact_cpu.inc](src/rff2/mb/verified_quality/compact_cpu.inc) — v202 | Local CPU adaptation of v201 `segments.hpp`, `disks.hpp`, `upper.hpp`, and the compact derivative design in `gpu_segments.inc`. It uses 192-bit MPFR centers with explicit outward errors, retains the full-P orbit, and retries full-P derivatives if compact verification is inconclusive. The creation notice identifies this independent application-code addition; it is not copied MPFR/GMP library source. Source/evidence: `Debug/source/v202/RESULT.md` and the included implementation. |
| [quality.hpp](src/rff2/mb/verified_quality/quality.hpp), [inverse_samples.hpp](src/rff2/mb/verified_quality/inverse_samples.hpp), [interval_solver.hpp](src/rff2/mb/verified_quality/interval_solver.hpp), [disks.hpp](src/rff2/mb/verified_quality/disks.hpp), [enclosure.hpp](src/rff2/mb/verified_quality/enclosure.hpp), [upper.hpp](src/rff2/mb/verified_quality/upper.hpp) | Retained local verification lineage, including the v81 family recorded during v112 integration and later disk/segment work already present at dd8f1c3. Some alternative solvers remain in redistributed source but are not the active proposal path. Their presence is not a claim that inverse interpolation or an interval-constraint solver now locates the nucleus. |


## Build references

The 2026-09-12 numeric update uses Merutilm's
[RFF-2.0 commit c305a990](https://github.com/Merutilm/RFF-2.0/commit/c305a990e71db016217afd3325eecceaa16f920c),
including the explicitly attributed GMP-derived multiply/square implementation
in [fixed_point_decimal.hpp](src/rff2/calc/fixed_point_decimal.hpp).
[fixed_point_complex.hpp](src/rff2/calc/fixed_point_complex.hpp) retains the local
caller API. Local additions handle tiny products and double exponent limits.

[FractalShark by Matthew Renzelmann](https://github.com/mattsaccount364/FractalShark)
and its reference-orbit/GPU notes were consulted for the later v214 orbit-reuse,
v215 NTT-fusion and persistent GPU experiments. Specific sections, source links
and the experimental scope are recorded in
[About the Algorithm, section 10.1](About_the_Algorithm.md#101-actual-implementation-and-early-project-sources).
They are separate from the numeric update and from the original v197-v202 work.

[BUILD_LOCATOR_SHADERS.md](BUILD_LOCATOR_SHADERS.md) records shader generation and compilation. The generator combines `ntt.comp`, `orbit_phases.inc` and `optimized_phases.inc`; CMake uses the host sources and eight compiled shader variants. [NFD_SOURCE_PROVENANCE.md](NFD_SOURCE_PROVENANCE.md) records the preserved Native File Dialog source revision.

License texts and component notices are maintained in [LICENSE](LICENSE), [NOTICE](NOTICE) and [THIRD-PARTY-LICENSE](THIRD-PARTY-LICENSE).
