<!-- Created by GPT-6 on 2026-09-12 -->

# Preserved Native File Dialog source

`extern/nativefiledialog-extended/` was copied from the actual local `build/_deps/nfd-src` working tree before that generated build directory was archived on 2026-09-12. Its recorded Git revision is `7bbbd9fe6b1d1549b41df138f614d1a44df9ba08`. The working tree reported no modifications before copying. The 192 retained source/document files were hash-compared with their original counterparts. Root Git metadata remains with the archived original; one nested submodule `.git` pointer copied initially was subsequently moved to `_local_only/nfd_git_metadata/`. No Git metadata is required to configure the preserved source.

The source reports Native File Dialog Extended version 1.3.0. Its original [LICENSE](extern/nativefiledialog-extended/LICENSE) and per-file notices are preserved unchanged. This local evidence records the source used by the prior v202 build; it is not a new online authentication of the commit or a rebuild of the library binary.

[CMakeLists.txt](CMakeLists.txt) now prefers this preserved source when present, using the same FetchContent dependency and target. The previous remote fallback remains available if the local source is absent. This selection change does not alter numerical code or the installed executable. A fresh CMake configuration is required after moving old build directories.

Copy hashes and the archived original path are retained in the repository's `_local_only/MOVE_MANIFEST.json`; that initial record includes the subsequently moved metadata pointer. `_local_only/NFD_METADATA_MOVE.json` records the pointer's final location and hash. Include this dependency source and its applicable notices when preparing the corresponding source package; do not treat it as disposable generated build output.
