<!-- Created by GPT-6 on 2026-09-12 -->
<!-- Modified by GPT-6 on 2026-09-12 -->

# Building the Locate Minibrot GPU shaders

These commands document the shader portion of the retained `Debug/source/v199/build_orbit.ps1` recipe. They were not executed during the license audit. The historical recipe names Vulkan SDK **1.4.335.0** and supplies no additional `glslc` optimization or target flags. This is a recipe record, not a fresh byte-for-byte rebuild claim.

Run from the directory containing this file in PowerShell. Use Python 3 on `PATH`; adjust the SDK path if necessary. The generator writes eight GLSL variants beside its source. The compiler writes the corresponding SPIR-V files there.

```powershell
$ErrorActionPreference = 'Stop'
$locatorShaderDir = Join-Path (Get-Location) 'src/rff2/mb'
$locatorGlslc = 'C:/VulkanSDK/1.4.335.0/Bin/glslc.exe'
python (Join-Path $locatorShaderDir 'make_orbit_shader.py')
if ($LASTEXITCODE -ne 0) { throw 'Orbit shader generation failed' }
foreach ($variant in 0..7) {
    & $locatorGlslc (Join-Path $locatorShaderDir "orbit_$variant.comp") -o (Join-Path $locatorShaderDir "orbit_$variant.spv")
    if ($LASTEXITCODE -ne 0) { throw "Orbit shader compilation failed: $variant" }
}
```

Retain all four generator inputs: `make_orbit_shader.py`, `ntt.comp`, `orbit_phases.inc`, and `optimized_phases.inc`. The output is `orbit_0.comp` through `orbit_7.comp`, then `orbit_0.spv` through `orbit_7.spv`. The generated GLSL retains its input notices.

After generating the SPIR-V files, configure and build the application using its Release toolchain and [CMakeLists.txt](CMakeLists.txt). CMake compiles `ntt_gpu.cpp` and `orbit_gpu.cpp`, and copies these eight SPIR-V files into `bin/` during configuration. It does not itself run the shader generator or compiler. Include the host sources, their included files, all generator inputs, and these instructions in the corresponding source package. Component notices are in [LICENSE](LICENSE).

The v199 recipe also compiles experimental test executables; those commands are not required for producing the application's shader artifacts and are not repeated here.
