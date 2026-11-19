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
- From LB (2026-11-19): a pan keeps the view's middle, not the pointer, so content under the pointer follows exactly only
  while no chip or gap enters the view (pass the pointer through `view.panPx`); a ring chip shows its label twice (the
  outside label and the pill); a block's label inline on its first text line (EC-95); `rows.ts`'s private `elide` onto
  `render/text.ts`'s.
- Pairs, ranges and maps on an alignment's columns (the author, 2026-11-19: "Yes, slice 4 (0.5.0)"; found by DT):
  `pairs.a`/`b`, `ranges.target` and `map.to` may name an `alignment` (its columns system), so a family's `SS_cons`
  arcs draw over the alignment and a structure's chain maps onto its closest seed row; in the fourth slice's 0.5.0
  bump with the molecules design. Until then the tRNA covariance view shows the seed's rows and logo, and the seed row
  as a sequence with its `SS_cons` pairs.
- AY's leftovers: tests for the rows and pairs layers' `items` (a row by its label, a pair by its arcs and dots); a
  dense figure's accessible layer about doubles the SVG (2,000 items in 388 kB of 732 kB): a coarser naming (one item
  a run of equal features) if a host minds.
- OP's leftovers: verify the bundled Liberation 2.1.5 files against the upstream release archive (they came from a
  LibreOffice runtime); check `image-rendering="optimizeSpeed"` in browsers (nearest-neighbour as resvg).
- DC's leftovers: `@molstar/io` as a dev dependency of seqx-examples so the record scripts read mmCIF with
  `parseCifText` (and the recorders out of the published `src/`); a ligand chip coloured on its own (one letters
  layer draws every block's pill).
- NB's leftovers: a text section's numbers in a source's numbering (inline only today); the `offset` and `segments`
  methods; a source's labels keyed other than by position; a records runtime converting labels of a declared
  numbering (files.md, Gaps).
- CY's leftovers: the label skip's time gain measured (the heavy frame-time lane has the scenario); a viewer test of a
  row request through `AddressModel`; `seq.resolve.length.unknown` before open quieter in the viewer; the 4 px
  arrowhead and the narrow-label constants in the config; Shift+drag across the origin (DP).
- DA's leftovers: the ANARCI oracle recorded by hand (the design's two rules checked); EC-29's indel carets (a map
  carrying insertions and deletions between two sequences, with 0.5.0's maps); short CDR labels at the fit view (H1
  by Kabat, L1 and L2 by IMGT); IMGT 112.1 never on the ruler at the closest zoom (three insertion labels a residue
  apart); the against view's gaps faint at the fit.
- QT2's leftovers: an `itWithin` option putting its unbounded half in medium (the plasmid pan at scale is written by
  hand, `@molstar/test-harness` shared: the author's call); core, viewer and storybook within 0.05 s of the 2 s guide.
- HC's leftovers: rulers naming their chain kept labelled (Complex, Antibody against); an empty track label still
  drawn as empty text (core `tracks.ts`); a cursor card replaced when a zoom re-hits the hover; the size script's
  `--files` crash.
- SF's leftover: look at a disabled toolbar button in `@molstar/ui`'s Storybook and in Prism, light and dark (the
  dimming computed at 1.7–1.9:1, not seen).
