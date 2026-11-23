# Drift: the keys for the next feature

**Where**: the viewer's keymap (`viewer/src/config.ts`, `view.keys`), CL3, item 4 (answer E9). **Found by**: CL3's
item 0 check, 2026-11-24.

**The design says**: answer E9 (`drafts/novice-questions.md`): `]` and `[` move the cursor to the next and previous
feature on the focused track, "Check that `]` and `[` are unbound"; the alternative Alt+Right and Alt+Left.

**The code does**: `]` and `[` are bound to the rows' zoom (`rows.zoom-in`, `rows.zoom-out`, EC-155; viewer spec
Bindings and Configuration). Alt+Right and Alt+Left are the view's forward and back.

**Free**: `.` and `,` (the ring already uses them for the next and previous partner of the focused stem), `n` and
`Shift+N`.

**Question for the author**: which keys walk a track's features?

**Recommended**: `.` and `,`: next and previous in the focused part, whatever its items are (a track's features, the
ring's partners), one meaning across the viewer; the rows keep `]` and `[`. Alternative: `n` and `Shift+N`.

**Decided by the author**, 2026-11-24: ". and , as recommended". `.` and `,` step to the next and previous item in the
focused part (a track's features, the ring's partners); the rows keep `]` and `[`. To CL3's item 4; the note is
deleted when CL3 merges with the answer in its design and spec.
