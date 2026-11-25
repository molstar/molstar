# Lane HP: drafting the fourth slice's track header parts design

**Task** (the author, 2026-11-19: "Fourth slice"; plan.md, Out of scope, "Track header parts"): draft the design by
[design-drafts.md](design-drafts.md) (read it all), on this topic: header actions as registry entries (icon, label,
when, command) drawn in the header's fixed slots, hit and announced; the built-ins (hide, color, filter) on the same
entries; a host adding or removing its own; a hidden-tracks chip in the section header that shows a hidden track
again; the headers' geometry as an API for a host's own widget later. Today: `track.hide`/`track.show`,
`layer.hide`/`layer.show` and `layer.filter` as commands, the "n filtered" chip the header's only control
([extensibility.md](../../../design/extensibility.md): "internal until an outside user"). Also read
`followups.md`'s notes on track header parts and legends' library surface, [ux-foundations.md](../../../design/ux-foundations.md)
(9. Chrome, the matrix's "a track or layer" row), [themes.md](../../../design/themes.md) (a track's colors, "color by"),
the structure syncing draft beside it (`drafts/structure-syncing.md`: its "color by this track" action), and the
precedents (the current Seq\*, IGV, JBrowse, Nightingale's headers).

**Output**: one file, `.scaffolds/seq/execution/third-slice/drafts/track-header-parts.md` (a fourth-slice draft). No other file,
no code.

**Model**: strong (opus). **Where**: the main checkout, read only except that one file. CF and MV work in worktrees
beside you.

**Report**: as design-drafts.md says (under 300 words).
