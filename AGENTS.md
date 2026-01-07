# Project rules for Codex (MATLAB)

- Do not rename any functions, files, structs, or fields used by the pipeline.
- Preserve all public function signatures exactly.
- Make minimal, surgical changes only where needed for correctness/stability.
- All new comments must be in English.
- Prefer adding small helper functions over large refactors.
- Always add/extend a fast smoke test when fixing runtime/selection bugs.
