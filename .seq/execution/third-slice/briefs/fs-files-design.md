# Lane FS, part 1: the files-as-sources design (a design only)

Follow `briefs/design-drafts.md` (Read first, What to write, Rules, Report) for one topic: **files as sources, and
the transforms a file's data needs**. Draft `.seq/execution/third-slice/drafts/files.md`. No code; no other file
edited. Model: strong (opus). Base: `packages-seq` at the commit this brief is in.

## The topic, in the author's words (2026-11-19)

- "we should have storybook example that load data from files too instead of just inlining the values. This includes
  the annotation data with descriptions of the transforms that need to be applied. Similar to how the base
  packages/sequences does that"
- "Also json sources so we can import uniprot stuff directly? And possibly other "well formatted" json?"
- "since we support json, it should support msgpack too out of the box, and both the SoA and AoS layouts for data"
- "Maybe the msgpack stuff could be moved to core/data/??? And reexported in io?" (decided: `base-library.md`'s
  MessagePack row: the encoder/decoder to `@molstar/core/data/msgpack`, the tree's codec to
  `@molstar/core/repr/codecs/msgpack`, `@molstar/io` re-exporting both).

## What it must cover

- **The tables text format** the sources design lists but no runtime reads yet: TSV and CSV (quoting, headers,
  comments, encodings); weigh GFF3 and BED (columns fixed by the format, 0- or 1-based by the format) as presets of the
  same reader, or not this slice.
- **JSON and MessagePack**: records found by a path (UniProt's REST entry JSON, `features[]`, `location.start.value`,
  `location.end.value`, `type`, `description`; the EBI Proteins API; any array of records); a UniProt entry's sequence
  and its features from one file; **both layouts**, an array of records (AoS) and an object of column arrays (SoA),
  read the same, typed arrays kept as columns without a copy.
- **The transforms**, written in the spec beside the file's source (never code): which columns or paths are the
  target, start and end; a 1- or 0-based start; an inclusive or exclusive end; the numbering positions are in; renamed
  and typed columns; a filter (the query sublanguage's, as `layer.data`); several annotations from one file (one per
  feature type). Read first: the base's `fields` remap, `from_source`, `uri` with `format` (`packages/sequences/_spec/
  format.md`, Annotation tables), and this codebase's sources, protocol and data-store designs (`.seq/design/sources.md`,
  `protocol.md`, `data-store.md`), and `packages/seq/core/src/data/` as built.
- **Where a transform runs**: in the runtime, once, producing the protocol's tables (the sources design: "formats
  converted at the edge"); its cache key; its issues (a column missing, a value not a number, a position past the end).
- **Examples and stories**: every dataset's fixtures loaded from files (FASTA, Stockholm, Newick, tables, JSON, the
  `seq` format) in the browser through the assets port and in Node through `node-assets.ts`; fixtures recorded by hand
  (a UniProt entry for P04637 saved, CC BY 4.0); tests never reach the network; fetching from services stays the
  showcases' registry entries (`ports.md`).
- **Minimal surface**: what is a spec param, a runtime option, a registry entry (`ports.md`'s rule); nothing configurable
  without a second use; host neutrality.

## Constraints

The language stays 0.4.0 this slice; if the transforms need a language change, say what and whether it waits for the
fourth slice's bump, or can be a source's params within 0.4.0. Name the `@molstar/core` additions (the MessagePack move
above; anything else is a gap).

## Report

As `design-drafts.md`'s Report, plus: the build lane's cut (one lane or two), and what the dataset lanes then do
differently.
