# Drift: what `packages`, `examples-and-stories` and `examples` say of the examples package and the folders

**What**:

1. **The examples package's dependencies.** `packages.md` (May import) and `examples-and-stories.md` say it depends on
   `@molstar/seqx-core` (and `@molstar/core`) only. It also depends on `@molstar/seqx-ext-rondo` (the genome ring's
   views) and has `@molstar/io` as a devDependency (the record scripts; the design's Tensions say so already).
2. **A view's build type.** `examples-and-stories.md` says `build(): Spec`; the modules return the builder's tree
   (`Node`, from `@molstar/core/repr/tree`).
3. **The licences list.** `examples.md`'s Configuration lists five ids (with `LicenseRef-anarci-numbering`); the code's
   `AllowedLicences` (`src/files/record-trim.ts`) has `CC0-1.0`, `CC-BY-4.0`, `LicenseRef-NCBI-public-data`,
   `LicenseRef-EMBL-EBI-terms`, `LicenseRef-seqx-generated` and `LicenseRef-UCSC-freely-usable`. The cytobands' licence,
   asked first in `ring-extension.md`, was settled as UCSC's "freely usable for any purpose".
4. **The folder graph is no test.** `packages.md` says each folder of a package is an entry, a module across folders
   names the entry, and the folder graph (no cycles) is `_test/layers.test.ts`. No such test exists in
   `packages/seq/core`; the workspace's import-graph check (`scripts/workspace/import-graph.mjs`) is about catalog and
   layer rules, not the package's folders.
5. **`packages/seq/AGENTS.md` does not list the examples package** (its table has core, viewer and storybook).

**Where**: `packages/seq/examples/package.json`, `src/views.ts`, `src/files/record-trim.ts`; `packages/seq/AGENTS.md`.

**Options**: 1 to 3: amend the designs (they are catalogs of intent; the examples spec now holds the rules as built).
4: (a) write the folder-graph test (a quick test over `src/*/` imports: each folder an entry, no cycles); (b) drop the
sentence. 5: add a row to the package guide (it is not a spec; the orchestrator's or the author's edit).

**Recommendation**: 1 to 3: amend. 4: (a), since the rule is cheap to hold and the spec's "an entry of its own" rules
(core's) depend on folders staying entries. 5: add the row.

**Note (fix round)**: items 1 and 2 are rows "design-only: not built as designed" in `packages` and
`examples-and-stories`; `examples › Dependencies` states what is built.
