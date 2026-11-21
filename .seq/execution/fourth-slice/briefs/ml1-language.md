# Lane ML1: the language: molecules, complexes, units, modified, compounds

**Task**: build the language half of [molecules](../../../design/molecules.md) (Decision: A position's chemical component, Molecules that are not linear polymers, Complexes and units, Links, Atoms to residues, On an axis, A section of a complex; Consequences, `seqx-core` `lang/`): the kinds `molecule`, `complex`, `unit`, `modified`, `compounds` with their params and checks, the builder's methods, issues `seq.lang.molecule.*`, `seq.lang.modified.*`, `seq.lang.unit.*`, `seq.lang.link.*`; the `seq` format's molecule, compound and complex tables; `source.ref` and `values.field` ([seq-data-ports](../../../design/seq-data-ports.md), Decided, gaps 2 and 4: a `source.ref` beside a `uri` or alone, alone an info issue, EC-185); `SeqLanguageVersion` 0.5.0 opened, you own the bump for the slice (print, patch, diff and the JSON checks cover every new field). In `@molstar/core`: the standard component per letter and alphabet (base-library.md's new row; `@molstar/core` public: state the surface in the first commit and stop for the author's acceptance). No resolve, render or runtime work (ML2, MR, MC).

**Output**: language, builder, format tables and core's constants with tests (round trips, checks, property tests on printing); `_spec` Language and Data; the Configuration rows of molecules.md that are spec params; the migration note for 0.5.0.

**Model**: strong (opus).

**Where**: A worktree on a branch `seq/ml1-lane` from the base the orchestrator names (after the third slice's close). Hot files (plan.md) are the orchestrator's: list the changes you need in the report. Every modified file's header updated; a quick-lane test beside the code, failing on the parent for a fix; the look in a visible pane, light and dark, before the report.

**Report**: Under 300 words: commits; checks run and their result; tests shown failing on the parent; the public surface added (names, for the author's acceptance where the brief says so); drift notes; what you saw; the `.v6/changelogs/` entry's text.
