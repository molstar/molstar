# Followups (draft; filed in packages/seq/TODO.md or core's TODO at finalize)

From lane A:
- `repr`: `checkRefs` (the resolve pass), `composeLanguage`, `patch` and `diff`, `descriptors`, issue dedupe (EC-86),
  catalog docs generation, a version-migration hook, `fallback` handling, a read-side `ParamValue` with typed columns
  only.
- The MessagePack codec in io runs `checkCodec`.
- A "Repr" section in `packages/core/TODO.md`.

From lane A's review (follow-ups, not blockers): union normalization re-validating slot subtrees (O(depth × size));
spelling checks of kinds and ids on trees built in code.

From lane B: `index.ts` re-exports `SeqLanguageVersion` as `LanguageVersion` (one definition); the core `_spec`'s language
section; a strand column check; codon exons against the protein's length (stop codons); a section's `range` revisited
with Locations.

From lane B's review: `repr` needs a version-aware validate (an unknown kind in a newer minor version is a warning,
format.md); the theme sublanguage slot for `layer.color`; the `refs` and `snapshots` levels in `SeqLanguage`.

From lane C: the `rect`, `row` and `chain` forms; lane E's static model implements `QueryModel` (the fixture tests become
conformance); per-key element access with access plans; EC-8 across a codon step (lane E); the `_spec` "Queries" section
(lane C's settled semantics: `not`'s universe, rect complement refused, zero-length overlap, snapping by label order,
missing ends as the system's ends).
