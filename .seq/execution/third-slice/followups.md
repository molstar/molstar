# Followups: the third slice

- `@molstar/ui`'s toolbar: a disabled button's color (`--msui-text-subtle`) is barely distinct from an enabled one's
  (`--msui-text-muted`); Seq*'s `mountToolbar` dims with opacity 0.4 meanwhile (C0). A fix in `packages/ui` would let
  that rule go.
- A toolbar button for `selection.bookmark` (C0).
- `@molstar/ui`'s `ToolbarItem.title` is both the accessible name and the tooltip, so the name carries the key
  (`Back (⌥←)`) (C0).
- `render/_test/frame-time.test.ts`'s max bound flaked once at the gate (18.8 ms against 16), passing on the rerun: the
  TODO item (best of three or a percentile).
