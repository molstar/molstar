# Building blocks

## Provisional

The whole note is the first proposal of the root nouns (§name-the-roots), drawn from the previous implementation's
concepts and meso\*'s language lessons. It is settled by the author's review and by the first examples running end to
end; until then names may change, and no lane builds on a name not yet used by a running example.

## The nouns

**Spec.** The root: `{ kind: "spec", params: { version, title? }, children: [data, views, refs?, snapshots?] }`. Data
is defined once; views refer to it; ids are local and every node, location, element and rect has an address
([addressing.md](addressing.md)).

**Coordinate system.** Where positions live. A system has a length (or is open, for a streamed genome), may be cyclic,
and has an alphabet when it carries letters. Positions are indices, 0-based, ranges half-open; labels (residue
numbers, insertion codes, chromosome names) are converted at the edge and printed back at the edge. A **structure**
system is a chain instance of a structure (an entry, a model, an assembly operator, a chain): its positions are the
chain's residues in the structure's own order, each keyed by the structure's numbering (its sequence index where it
has one, the author's number and insertion code) and flagged observed or not. It names the structure abstractly; a
structure runtime (Mol\* in practice) binds it to atoms, so the core never reads a structure file.

**Sequence.** A coordinate system with letters: a protein, a nucleic acid, a chromosome. Shared by its **units** (a
chain, a copy, a haplotype); units are grouped by a **complex**. The copy rule: a location on one unit never reaches
another unit through the sequence they share.

**Map.** A directed relation between two systems. Three questions, answered separately:

1. **What it means: its kind.**

   | Kind | Meaning | Compact form |
   |---|---|---|
   | `identity` | position k is position k | none |
   | `offset` | k is k + d | one number |
   | `numbering` | positions to labels by a scheme (Kabat, IMGT, …) or a label list | a scheme name, or labels |
   | `segments` | runs: from, length, to, direction | the runs |
   | `pairwise` | an alignment of two sequences | a CIGAR string |
   | `alignment` | a row onto its alignment's columns | the gapped row itself |
   | `codon` | residues onto nucleotide triplets | exon ranges and a strand |
   | `structure` | a sequence onto a structure's chain residues (SIFTS in abstract form) | segments plus residue keys, with observed flags and per-segment provenance |

2. **Where it comes from: its sourcing.**
   - **inline**: written in the spec, in its kind's compact form;
   - **a source**: a table or file read through the source protocol, by region when large (a chain file, a PAF, a
     SIFTS table, a GFF's CDS lines);
   - **computed**: produced by a runtime the application binds when it composes the viewer (align a chain to its
     UniProt sequence, look up SIFTS for an entry, translate a gene model into its codon map). A computed map keeps its
     kind: a codon map computed from a gene model is a `codon` map. The spec names the method and its inputs; without
     the runtime the map is unavailable, with an issue.

3. **How it is held: its encoding.** A resolved map answers one interface (the image of a range, forward and back, as
   runs; the boundary rules; the status of each element), and is held in the cheapest encoding that answers it:
   - **implicit**: arithmetic, nothing stored (`identity`, `offset`, a codon map within one exon);
   - **runs**: segments (a CIGAR, exons, an alignment row's gap runs), often a handful of rows;
   - **expanded**: one entry per position, only where nothing shorter holds (arbitrary labels, a structure's residue
     keys with gaps and insertion codes), and only for the region asked.
   An encoding is never expanded to answer a question a shorter one can: composing a path composes encodings (two
   offsets give an offset; runs through runs give runs), and a genome-scale map is held by region.

Direction (+1 or −1) is part of every run, so an inversion, a reverse-strand gene and a protein on its codons are one
mechanism (the JBrowse reply, point 3). A **path** composes maps (bounded length); **projection** carries a location
along the shortest unambiguous path, and reports when there is none. Structure maps make a structure one more system on
those paths: a variant on the gene reaches a residue in the structure through codon and structure maps, and the copy
rule keeps each chain instance apart.

**Alignment.** Rows (each a unit or a sequence with its map onto the columns) over a columns system.

**Annotation.** Data on a system, in one of six shapes: `ranges` (intervals with fields), `values` (a number per
position or bin), `pairs` (two locations, possibly on two systems: contacts, base pairs, synteny), `matrix` (a value
per row and position), `grid` (a value per pair of positions), `rows` (keyed records). An annotation reads from a
**source** ([sources.md](sources.md)).

**Location.** A set of ranges on one system, or of rects over two axes (positions × positions, rows × columns, states ×
positions): what a query or an address returns, what hover and selection hold, what focus names.

**Ref.** A user-defined name for an address, in the `refs` container: a bookmark, a named region, a set; usable wherever
an address is.

**Query.** A small composable language whose result is a location: generators (a ref, all, items of an annotation,
letters matching a motif), modifiers (range by label, expand, onto another system, group), filters, combinators. Its
canonical form is a tree like the spec's, encoded the same ways; it prints as S-expressions.

**View.** How data is shown: a **section** is one axis (a system, or a piece of one: a range, reversed, with collapsed
gaps) with its **tracks**; a track stacks **layers**. A view lays its sections out (`tracks` or `text`) and may name
the **focus** it opens on.

**Layer.** One way of drawing one annotation (or the letters) on a track: `letters`, `ranges`, `values`, `pairs`,
`matrix`, `logo`, `ruler`, … Each layer type declares its params and its **fallback** (the simpler layer a reader that
lacks it draws instead).

**Theme.** A channel (color, size) and how its value is computed: `uniform`, by a field's values (a palette, binned or
continuous, a domain open to the data in view), by rules (thresholds), by a named scheme. Data, never code.

**Snapshot.** A named patch over the spec: a story's step, a bookmark.

**Issue.** `{ code, severity, path, message, hint? }` with a stable namespaced code, from every pass (validate,
resolve, load, draw); reported, never thrown, unless the spec cannot be read at all.

## How they fit

The levels these nouns live in, and the IRs a spec compiles to, are [composition.md](composition.md).
