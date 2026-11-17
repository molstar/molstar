# Followups (draft; filed in packages/seq/TODO.md or core's TODO at finalize)

From lane A:
- `repr`: `checkRefs` (the resolve pass), `composeLanguage`, `patch` and `diff`, `descriptors`, issue dedupe (EC-86),
  catalog docs generation, a version-migration hook, `fallback` handling, a read-side `ParamValue` with typed columns
  only.
- The MessagePack codec in io runs `checkCodec`.
- A "Repr" section in `packages/core/TODO.md`.

From lane A's review (follow-ups, not blockers): union normalization re-validating slot subtrees (O(depth × size));
spelling checks of kinds and ids on trees built in code.
