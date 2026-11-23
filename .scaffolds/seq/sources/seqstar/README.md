# The previous implementation (packages/sequences)

Seq\* as built through 2026-11: `@molstar/seq-core` (format, model, query, builder, streams), `@molstar/seq-viewer`
(headless context, Canvas2D and SVG backends, React and DOM UIs), a Nightingale adapter and a Storybook. About 103k lines
of source and 119k of tests. Its `_spec/` (README, model, format, seqql, viewer, interaction, sources) is the best
statement of the concepts; read it for meaning, not for structure. It will be deleted once `packages/seq` replaces it,
so `packages/seq` never names it.

## What we keep (the concepts)

- **Positions are indices**: 0-based, half-open ranges; labels (`"27A"`, insertion codes, numbering schemes) are
  converted once, at load. A ref is identity.
- **Coordinate systems and maps**: every map is a segment table (position, length, target, direction), with a codon
  step of span 3; maps may be partial; a path composes mappings (bounded length, at most one codon step); projection
  carries a location from one system to another; the copy rule (chain A never reaches chain B through a shared
  molecule).
- **Molecules, units and copies**: a molecule (sequence) is shared by its units (chains, copies); a complex groups units.
- **Annotation geometries**: ranges, values, pairs, matrix, grid, rows.
- **Data defined once, views over it**: `data` holds sources, sequences, alignments, mappings, annotations; `view`
  holds sections, tracks and layers that refer to data by ref.
- **Layer types with fallbacks**: a reader that does not know a layer type draws its fallback (down to letters, ranges
  or nothing) and says so.
- **Themes as data**: a channel (color, size, glyph, line style) takes a theme by name with params.
- **Patches with inverses, snapshots as patches over the base.**
- **A headless context** that lays out, hit-tests and draws through an abstract backend (Canvas2D, SVG export).
- **Hover, selection and highlight as locations**, projected onto every section that shows the same positions.
- **Behaviors as replaceable policy**, commands for everything the viewer does.
- **Issues with stable codes** (`{code, severity, path, message}`), from every phase.
- **SeqQL**: a composable query language (generators, modifiers, filters, combinators) whose result is a set of
  ranges.

## What we leave (accidental complexity)

- **Seven test tiers** (quick, perf with byte-exact allocation baselines, heavy, browser, heap, browser-heap,
  heap-deep): about 190 guard files, and most of its AGENTS.md about V8 tier-up flags. Here: the repository's lanes.
- **Zero allocation everywhere** (struct-of-arrays tables, out-parameters, module scratch, enums hoisted to module
  constants): essential only in the draw and hover loops; elsewhere plain code.
- **Process scaffolding in the code**: phase ids in comments, about 900 spec citations, a spec checker for citations,
  headers and guard templates, line budgets, generated reference tables written by snapshot tests.
- **Surface area**: two default UIs (React and DOM) plus Nightingale; 21 behaviors with exact binding params; about 25
  managers on the context; a renderer interface with more than 25 optional members; document layering
  (base, snapshot, local, preview) with masked patch routing.
- **Streaming as a host-specific tile protocol** with opaque stream ids, and rendering semantics (value domains,
  palettes, sequence lengths) that depend on a whole-data manifest. The JBrowse reply shows why this fails
  portability ([jbrowse-reply](../jbrowse-reply/README.md)).
- **rxjs** (`SeqState`, `Subject` streams, `Subscription`) and its own primitives beside it (`Live`, `SeqEvent`,
  `Fanout`, a scheduler port, a dependency graph). Here: `@molstar/core/reactive`.

## Examples it had

22 cases (p53 text and tracks, antibody numbering, protein family, RNA base pairs, DNA construct, imatinib in ABL,
barnase–barstar, AFDB complex, PF00042 alignment, story drafts, designed dimer, UniProt features, protein over coding
DNA, indels in wild-type numbering, plasmid ring, genotype matrix, queries and bookmarks, streamed chromosome, protein
on its genome through codons, reverse-strand gene, family triage) and topic stories (alignments, axes and strands, dot
plots, navigators, themes, hover, patches, SeqQL, search, streaming, performance, text view, file formats). The new
examples are fewer and each shows one building block ([design/examples.md](../../design/examples.md)).

## Reconciliation

- Its "SoA everywhere" contradicts the reading of §performance as designed-in where it matters; kept for hot loops
  only.
- Its document layering is transient state (§committed-vs-transient) built as a stack of documents; here a preview is
  a transient overlay over one committed spec.
- Its spec-to-code citation machinery is replaced by the root's rule: code comments cite the `_spec` section they
  follow, nothing more.
