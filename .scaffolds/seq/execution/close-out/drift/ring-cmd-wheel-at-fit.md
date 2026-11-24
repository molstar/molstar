# Drift: Cmd/Ctrl+wheel over the genome ring at the fit

**Level**: code-vs-design. **Outcome**: (c), decided by the author (2026-11-24, "Agree all"); to FX27. FX26 already keeps the hidden window still on a ring nobody follows. **Found
by**: CL12, 2026-11-24 (M5's gap 2; novice 4).

**What**: `layouts.md` (Rings): "A ring nobody follows has no window: its zoom and pan are idle, as a navigator's
(EC-110)"; the matrix row "a ring section": zoom "idle with no follower". `ring-extension.md`'s matrix, "a section: the
genome ring", zoom: "no: a ring shows its whole axis; opening a chromosome is `section.include`". Neither says whether
an idle zoom is taken (the page does not scroll, nothing changes) or given back to the page. The brief asks for never
"taken but nothing". The code does a third thing.

**Where**: `viewer/src/input.ts` `wheel` (the binding table's `zoom` taken whenever it resolves), `viewer/src/behaviors/
navigate.ts` `wheel` → `ViewModel.zoom`; `viewer/src/models/view.ts` `_zoomX` (no idle check for an unfollowed ring).

**Evidence** (headless, Features/Genome Ring: The Genome At The Fit, one ring section, no follower, `canMove` true):
- Cmd+wheel in over a stem: taken (`true`); the section's viewport goes from 0–3,031,042,417 to 568,192,054–
  2,993,025,988. The ring still draws the whole genome: an invisible window moved. Cmd+wheel out: back to the whole.
- A plain wheel at the fit: taken, nothing moves (FX24 kept it taken: a pan never chains to the page).
- So the novice's "taken but no change" is a hidden window narrowing, against "no window: idle".

**Options**:
(a) Idle and not taken: no viewport change, the wheel returns false, the page gets it. On Windows and Linux Ctrl+wheel
then zooms the browser page; on macOS a pinch (Ctrl+wheel) zooms the page too. Honest, but surprising.
(b) Idle and taken, said: no viewport change; the first such wheel of a gesture announces once, at the pointer or in the
live region, "The ring shows the whole genome: open a chromosome from its label's menu or the square" (no repeat while
the wheel keeps going).
(c) The ring's zoom is the cut: Cmd+wheel in over a chromosome opens it (`section.include` of the block under the
pointer, one patch, one undo step), out on an opened chromosome shows all; one cut a wheel gesture (deltas past a notch,
then nothing until the wheel rests). The designs already call opening a chromosome "the ring zoom" (`ring-extension.md`,
The ring zoom: a cut); this gives it the zoom gesture.

**Recommendation**: (c), with (b)'s words over the hole (no chromosome under the pointer); and in every case the hidden
window no longer moves on a ring nobody follows (`_zoomX` and `pan` idle there, as `layouts.md` says). (b) alone if the
author wants the wheel never to change which blocks are shown.
