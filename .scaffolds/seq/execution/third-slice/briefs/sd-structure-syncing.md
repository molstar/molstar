# Lane SD: drafting the fourth slice's structure syncing design

**Task** (the author, 2026-11-19: "Structure syncing", for the fourth slice; plan.md, "Structure integration in the
fourth slice is a fake 3D viewer"): draft the design by [design-drafts.md](design-drafts.md) (read it all: what to
read, the nine sections, the rules, the report), on this topic:

- A small storybook component stands in for a structure viewer: its chains and residues as a schematic, taking
  MolViewSpec-shaped states and sending hover, selection and focus as a structure viewer would. It shows every syncing
  path: hover and selection both ways, focus, "color by this track", a structure's chains mapped onto sequences
  (copies, author numbering, unobserved residues). The ports, maps and messages are the real ones; the Mol\*
  integration (Mol\* `Model` adapters, a live 3D view) waits for the new mol viewer that replaces Prism ("so i don't
  want to waste work").
- Designs to read and build on, not restate: [ports.md](../../../design/ports.md), [maps.md](../../../design/maps.md)
  (structure maps, the hull), [molecules.md](../../../design/molecules.md), [composition.md](../../../design/composition.md),
  [viewer-composition.md](../../../design/viewer-composition.md) ("Inside a Mol\* plugin"),
  [locations-and-queries.md](../../../design/locations-and-queries.md), [themes.md](../../../design/themes.md) (a
  track's colors), [examples.md](../../../design/examples.md) (p53's 2OCJ, Complex's 1FIN), `followups.md`'s
  "seq data ports" and "track header parts" notes (what they leave to this design), MolViewSpec
  (`packages/mvs/`: its state shape and its selectors), and how the current Seq\* syncs with Mol\* (precedent only).
- Host neutrality is the rule to hold: Seq\* imports nothing of a host; the fake viewer is a host and talks only
  through ports and messages a real one would use.

**Output**: one file, `.scaffolds/seq/execution/third-slice/drafts/structure-syncing.md` (a fourth-slice draft, moved to
`design/` when the author agrees). No other file, no code.

**Model**: strong (opus). **Where**: the main checkout, read only except that one file (no installs, builds or git
writes). CF and MV work in worktrees beside you.

**Report**: as design-drafts.md says (under 300 words).
