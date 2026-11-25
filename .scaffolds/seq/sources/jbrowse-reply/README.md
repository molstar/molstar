# The JBrowse reply

A reply from the JBrowse side to the Seq\* developers (the seqviewspec-reply exploration repository, dated
2026-11-16, against SeqViewSpec 1.0.0). Its question: could JBrowse open a SeqViewSpec document with its own adapters
and renderer? Its answer: yes, for most linear-genome content; synteny, dot plots, circular views and Hi-C with help
beyond the kernel. Each claim has a probe document and the validator's message. Friendly, empirical, framed as
questions. It carries eleven examples, a spike converter (document → JBrowse config) and a gallery.

## What it asks, ranked by how much it changes the design

1. **Self-describing sources.** A document can name only the host's tile protocol and opaque stream ids; once the
   JBrowse-private section is stripped, every example draws nothing. A document should say what the data is: its
   format (bigWig, bigBed, 2bit, VCF with a tabix index, FASTA with fai, GFF3, BAM, CRAM, PAF), its URI and index,
   the assembly's identifiers and aliases, and which file column is which field. *Taken, with the author's twist: one
   protocol, formats converted onto it by runtimes the application provides ([design/sources.md](../../design/sources.md)).*
2. **Fetch by region; nothing depends on a whole-data manifest.** `domain: null` should mean the data in view; color
   keys stable without a global vocabulary (keep the FNV-1a hash for unlisted values); sequence names and lengths come
   from the document or the file, not a manifest. *Taken.*
3. **Strand on every mapping and on pairs**: an inversion between genomes has no form today (only codon mappings take
   a strand). *Taken: direction is part of every map segment.*
4. **Layers across sections**: ribbons joining two sections, a dot plot with its own y section; streamable pairs,
   matrix and binned grid data. *Deferred; the coordinate model must not rule it out.*
5. **A published JSON Schema and wire contract**, so a reader needs no 3.5 MB core. *Deferred by the author (2026-11-17: a schema may bloat, and specs are mostly built with the builder); the
   description could derive one later. Was: generated from the
   language description.*
6. **Unknown kinds don't fail the document** (they become unavailable, with an issue); `diff(a, b)` returns the patch
   from a to b. *Taken.*
7. Smaller: `split` on values, a log scale, a points encoding, numeric filters (min, max), lane order on ranges,
   joining a layer to a keyed table, heights in px, multi-level gene hierarchies, mapping files (chain files) by
   reference, derived-ref characters for HLA names. *Deferred to the examples that need them.*

## What it says to keep

The kernel/extension split with fallbacks ("made even a quick loader practical"); the FNV-1a hash for unlisted values;
threshold `rules` themes, diverging palettes, shared domains, binned palettes; liftover through segment mappings; node
shape and refs; 0-based half-open indices; rows sharing an `id` as one item; `default_view`; `direction` on sections and
blocks; group and summary semantics; the per-type fallback table; the `complement` letter map.

## Reconciliation

- Points 1 and 2 are §bulk-data-not-nodes and ports (architecture): data is referenced, a runtime resolves it.
- Point 5 follows from languages as data ([meso](../meso/README.md)): a schema would be one more derivation, deferred.
- Point 6 is §parser-repairs-slips' honesty: report, never guess, never fail what can be skipped.
