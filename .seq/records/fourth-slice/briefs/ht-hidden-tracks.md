# Lane HT: hidden tracks are the application's

**Base**: `packages-seq` at the orchestrator's head after FX6 (`git reset --hard packages-seq` first). Branch `seq/ht`.
Model: fast. The fix lane's rules ([fx-fix-lane.md](fx-fix-lane.md)) on headers, tests and hot files hold.

**Task** (the author, 2026-11-21: "hiding a track moves the UI around … give this to control of the app layer and that
is then responsible to showing the hidden tracks and ability to restore them"): build track-header-parts.md's amended
rule ("Hidden tracks are the application's", the section header, chips, the entries table, Configuration, Tensions):

- the view draws nothing for a hidden track: no section-header line, no "n hidden" chip; `section.header: auto` only for
  a slotted section entry; core's built-in header states lose the section's hidden-tracks chip;
- the viewer exposes each section's hidden tracks (from the spec, in order: address, label, section) with a change
  event, and `track.show` takes one, several or all (one patch, one undo step);
- `seq:hidden-tracks` becomes a menu-only section entry ("Show hidden ▸": each hidden track, then "Show all");
- `mountHiddenTracks(viewer, el, options?)` in the viewer (beside `mountLegend`; `@molstar/ui`'s list), a host places it;
  the default viewer's Storybook composition mounts it in its toolbar (beside the back/forward and keys buttons), the
  minimal viewer none;
- the core and viewer `_spec`s and their Configuration tables; the stories (`Features/Track headers`, the story that
  showed "1 hidden") updated; their twins.

Tests: hiding a track leaves every other track's y unchanged above it and moves those below by exactly its height (no
line added); the list shows and restores one and all; the menu's "Show hidden ▸" by keyboard; a section with all
tracks hidden collapses and a jump into it says "on a hidden track". Look in light and dark on port 6014.

**Checks and report**: as [fx5-author-look.md](fx5-author-look.md) lists them; under 200 words, the public surface
added and removed, the `.v6/changelogs/` text.
