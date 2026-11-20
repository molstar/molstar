# Lane MK: structure maps: segments, priority, via

**Task**: build the map half of [structure-syncing](../../../design/structure-syncing.md) (Choosing between two maps of one pair of systems; Consequences `resolve/`, `lang/`) and [maps](../../../design/maps.md): `segments`, the `structure` kind with per-residue status, a map's integer `priority` (default 0; distinct shortest candidates ranked by their weakest map; a tie stays ambiguous; the carried search dropped when it changes) and a section's `via` (a ref or a list; rivals of the named maps left out for that use; `follow`, `direction: { along }`, projection of hover, selection and the focus hull), `StaticModel.path(from, to, via?)` keeping one entry per pair and `via` set, `seq.resolve.path.ambiguous`'s hint naming both; a numbering's `offset` and `segments` methods (NB's leftovers). Both params go through ML1's bump, added to its language section. EC-22 to EC-27 on implicit, run and expanded encodings with a ranked pair, a tie and a `via`; EC-24 amended as the design says (unless ranked or named).

**Output**: maps, paths and their tests (a priority never lengthens a path; `via` leaves other paths alone); `_spec` Resolve (Paths) and Language; maps.md's Design page text if it moved.

**Model**: strong (opus).

**Where**: A worktree on a branch `seq/mk-lane` from the base the orchestrator names (after the third slice's close). Hot files (plan.md) are the orchestrator's: list the changes you need in the report. Every modified file's header updated; a quick-lane test beside the code, failing on the parent for a fix; the look in a visible pane, light and dark, before the report.

**Report**: Under 300 words: commits; checks run and their result; tests shown failing on the parent; the public surface added (names, for the author's acceptance where the brief says so); drift notes; what you saw; the `.v6/changelogs/` entry's text.
