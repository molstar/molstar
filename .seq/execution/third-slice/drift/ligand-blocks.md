# Drift: a ligand block in tracks and rings

**What** (the author, 2026-11-19, Complex › Blocks › Row, light: "the ligand rendering really isn't good ... we need to
find a better way"): ATP, a one-position block after the polymers, gets one position's share of the row: at the fit view
a sliver at the right end, its label cut to "A"; zoomed in, "AT" at the edge and a grey cell in the sequence track. The
designs give a one-position block a form only in text (layouts, Blocks, "Text": chips that flow like words, EC-99);
tracks (`arrange: row` and `stack`) and rings have no rule, so it falls to the linear share.

**Where**: `.seq/design/layouts.md` (Blocks: Tracks, Rings); `render/layout/` (the row's block widths, the gaps' share,
EC-102); the Complex and Globins views (ATP in 1FIN; the hemes in 4HHB, EC-64).

**The base (`packages/sequences`)**: chips in its text view (the same rule we kept), and in tracks "a ligand as a pill
over its cells, its label fitted; narrower than its box, a filled cell at least 2 device px wide, unlabeled, so it never
vanishes". So in tracks it had the same problem, floored: at the fit view a 2 px unlabeled mark. Not better where it
matters.

**Options**:

- **(a) A chip in every layout.** A one-position block takes a fixed slot in screen px, its label's width (the text
  chip's rule: the component and `×N`, 5 px either side), like a gap and unlike a position: zoom does not scale it, so
  it reads at the fit view and never clips. Tracks: the slot after the polymers' blocks (the row) or a short row of
  chips under the stack; the block's tracks draw one pill a track holding something there (the letters' band its
  label, a feature band a mark), empty bands blank. Rings: a fixed-angle arc after the polymers' arcs. One rule for
  text, tracks and rings; hit, hover, selection, links read it as its one position (as the base). Cost: the row's x is
  no longer linear over blocks, as with gaps already (EC-102); `windowExtent` and the ruler skip the slot as they skip
  gaps.
- **(b) Ligands off the axis.** A chip rail under the section's labels (ATP ×2, heme ×4), not blocks; hover and select
  as positions; the axis is polymers only. Simpler layout, but a ligand stops being a block (the language's `block`,
  EC-101 change) and has no place in a track's columns.
- **(c) The atoms form.** `of` the per-atom system: ATP 31 positions wide. Already planned as a choice (Complex views);
  as a default, it spends width on letters (`X`, elements) that mean nothing to most readers.
- **(d) Where it binds, on top of (a).** The ligand's contacts drawn as a `pairs` layer from the chip to the residues it
  touches (K33, E81, L83 in CDK2), the builder's "contacts" ref highlighting them: the reader sees not only that ATP is
  there but where it sits. Already planned in the Complex and Globins views (a `pairs` layer of the derived contacts);
  (a) gives the arcs a place to land at every zoom.

**Recommendation**: (a) as the layouts' rule for every one-position block, with (d) as the Complex's and Globins'
default view of their ligands (the contacts layer PG builds). (c) stays a choice per block. Not (b): it changes the
language for a drawing problem.

**Scheduled**: brought to the author now (2026-11-19); once decided, layouts.md and its Design page updated, and built
by lane LB in stage 5, after SO and before DG and DC (both draw ligands: the hemes, ATP).
