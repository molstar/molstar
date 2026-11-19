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
  headers' geometry for a host's widget later (plan, Out of scope). Its story: "color by this track" on a fake 3D viewer
  (the Mol\* integration waits for the new mol viewer, decided 2026-11-19), a click building a MolViewSpec state ad hoc.
- Legends' library surface, for the fourth slice (the author, 2026-11-19): display, placement, hiding and collapsing,
  a compact form; editing from a legend is the application's, on the commands (plan, Out of scope).
- For the fourth slice (the author, 2026-11-19): seq data ports (representation ports) (share Mol*'s parsed representation,
  not only bytes); custom navigators (a navigator's picture as a registry entry, the window behavior public); event
  composition across parts is the application's (plan, Out of scope).
- Structure integration, fourth slice (decided 2026-11-19): a fake 3D viewer demonstrates the syncing (hover, selection
  and focus both ways, "color by this track", chain maps) on the real ports and messages; the Mol\* integration and
  `Model` adapters wait for the new mol viewer that replaces Prism (plan, Out of scope).
- Theme clause kinds as symbols (the author, 2026-11-19: "we could support optional aliases with symbols for = etc if
  needed later"): optional aliases (`=`, `<`, `<=`, `>`, `>=`) read as the clause words; only when a user asks.
- A minimal viewer's size (the author, 2026-11-19, on the size guides: "i guess that's fine. this is for the default
  composition of the viewer? i think a minimal one could be smaller"): measured in stage 7 (M3).
- Query forms the design names and QS found unbuilt (2026-11-19): `search` (a motif in the letters), `onto` (a result
  carried to another system as a query), `expand`, `group`; Features/Query gains their stories when built.
