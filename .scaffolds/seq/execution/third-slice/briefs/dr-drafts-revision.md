# Lane DR: the two fourth-slice drafts revised by the author's decisions

**Task**: revise `drafts/structure-syncing.md` and `drafts/track-header-parts.md` by their "Decided by the author"
sections (appended 2026-11-19), as [design-drafts.md](design-drafts.md) shapes a draft: the Decision, the matrix
cells, the Configuration table, Tensions, the self-check and Consequences changed to match; the gaps answered moved
out of "Gaps for the author". New for structure syncing: Seq\* declares no MolViewSpec expressions (the port speaks Seq\*'s
locations; the host's adapter translates; the fake viewer's adapter shows it); the `structure` map kind kept; and a
design for **choosing between two maps joining the same pair of systems** (today ambiguous, EC-22–27: e.g. SIFTS and an
alignment for one chain): read [maps.md](../../../design/maps.md) and `core/src/resolve/paths.ts`, propose how a section
or link names the map it goes through (or a map's priority), as one new gap with a recommendation. New for track
header parts: both the canvas-drawn slots and a host's DOM parts over the header geometry; the secondary click a
behavior param.

**Output**: those two files only. No code.

**Model**: strong (opus). **Where**: the main checkout, read only except those two files.

**Report**: under 300 words: what changed in each, the new gap with its recommendation.
