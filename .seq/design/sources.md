# Sources: one protocol, formats converted onto it

## Intent

A spec must say what its data is, so it opens anywhere: no opaque ids only one host can resolve, no meaning that
depends on a manifest of all the data. And the viewer must read one thing, not a parser per file format.

## Precedents

- The previous SeqViewSpec: a source named a URI, a service, a host file, a value, or an opaque `stream` id served by
  the host's tile protocol (`seqstar-tiles/1`); value domains, palettes and sequence lengths came from a whole-data
  manifest ([sources/seqstar](../sources/seqstar/README.md)).
- The JBrowse reply: once the host-private section is removed, every example draws nothing; a document should name the
  format, the URI and index, the assembly and the column-to-field map ([sources/jbrowse-reply](../sources/jbrowse-reply/README.md)).
- JBrowse itself: one adapter per format, all yielding features by region.

## Decision

The author, 2026-11-17: "i still want this to have a single protocol. but what the spec should have is the ability to
define a source in a format which an application provides a "runtime" for .. meaning the streaming protocol should be
standardized, and we convert other formats onto it".

- **A source names its data in the data's own terms**: a `format` (`fasta`, `bigwig`, `2bit`, `vcf`, …, or `seq` for
  data already in the protocol's shape), where it is (`uri`, plus `index` when the format has one; or an inline
  `value`), and what it holds when the format does not say (the assembly's identifiers, a column-to-field map). The
  spec is complete without host code.
- **Answers carry nodes as well as tables** (sequences, maps, annotations, nested sources), identified by (source,
  key) ([dynamic-nodes.md](dynamic-nodes.md)).
- **One protocol** is what the model and the viewer read: a source, opened, answers *what it covers* (its sequences
  and lengths, or the region it can serve) and *give me this region at this level of detail* as typed tables in the
  model's own shapes (sequences, ranges, values, pairs, matrix, grid, rows). Small inline data is read through the same
  protocol, so there is one path.
- **Runtimes bind to resolved sources**: a source resolves to an entry of the resolved model (its format, location, coverage), and the
  application's runtime for that format binds to the entry, never to the spec node ([composition.md](composition.md)).
- **A runtime per format**, the format's default adapter of the seq data ports ([seq-data-ports](seq-data-ports.md)):
  a host's adapter replaces it for the entries it holds in memory, one adapter a format; provided by the application: a registry entry the viewer's spec lists (§executors-by-composition),
  that turns a format into the protocol. The core ships the protocol, the `seq` format and the small text formats
  (FASTA, alignments, tables, Newick; the tables, JSON, MessagePack and GFF3 through one records runtime whose
  transform is params of the `source` node, [files.md](files.md)); binary and indexed genome formats are runtimes an
  application adds. The genome runtimes ([genome-formats](genome-formats.md), the fifth slice), each an opt-in entry
  of `seqx-core` like `runtimes/mmcif` and `runtimes/bigwig`: `runtimes/bed` (BED, bedGraph, BEDPE, chrom.sizes; built, BE),
  `runtimes/gff` (GFF3 and GTF as transcripts) and `runtimes/vcf`, reading BGZF and tabix from `@molstar/io`; the `seq`
  runtime serves sided pairs from either end and a map's runs. Designed only, no runtime (a host brings its own through
  `withAdapter`): bigBed, 2bit, FASTA by region with `.fai` or `.gzi`, PAF and PIF, chain files (genome-formats,
  Appendix). A source
  whose format has no runtime is unavailable, with an issue; the rest of the spec opens.
- **No whole-data manifest in the meaning**: a value domain left open means the data in view, computed when the view settles ([data-store.md](data-store.md)); a categorical value not
  listed in the theme takes a color by a stable hash of the value; lengths come from the source's coverage answer.

## The protocol's types

Settled in [protocol.md](protocol.md): runtimes, coverage, requests, answers and their table shapes, errors, the `seq`
format, the store's keys and request planning across maps.

## Rejected alternatives

- Opaque stream ids resolved by the host (the previous design): not portable.
- A list of registered protocols, one per format, read by the viewer (the reply's proposal as written): the viewer would
  grow a reader per format; here formats are converted at the edge, the author's choice.

## Consequences

- The model's annotation shapes are the protocol's answer shapes: one set of tables, no conversion in between.
- Conformance: every runtime passes one suite (coverage, region answers, empty and edge regions, cancellation), run
  against fixtures without the network.
- A runtime lands with its synthetic twin: a source of the same protocol whose answers are delayed on an injectable
  clock, which its stories and tests use ([testing.md](testing.md), How tests are written).
