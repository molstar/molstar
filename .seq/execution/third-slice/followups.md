# Followups: the third slice

- `@molstar/ui`'s toolbar: a disabled button's color (`--msui-text-subtle`) is barely distinct from an enabled one's
  (`--msui-text-muted`); Seq*'s `mountToolbar` dims with opacity 0.4 meanwhile (C0). A fix in `packages/ui` would let
  that rule go.
- A toolbar button for `selection.bookmark` (C0).
- `@molstar/ui`'s `ToolbarItem.title` is both the accessible name and the tooltip, so the name carries the key
  (`Back (⌥←)`) (C0).
- `render/_test/frame-time.test.ts`'s max bound flaked once at the gate (18.8 ms against 16), passing on the rerun: the
  TODO item (best of three or a percentile).
- Track header parts, decided for the fourth slice (the author, 2026-11-19): header actions as registry entries in
  the header's fixed slots (a host's button among them), a hidden-tracks chip that shows hidden tracks again, the
  headers' geometry for a host's widget later (plan, Out of scope). Its story: "color by this track" with the Mol*
  integration, a click building a MolViewSpec state ad hoc for the 3D view.
- Legends' library surface, for the fourth slice (the author, 2026-11-19): display, placement, hiding and collapsing,
  a compact form; editing from a legend is the application's, on the commands (plan, Out of scope).
