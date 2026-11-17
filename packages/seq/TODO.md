# TODO

What is not built. An item leaves in the commit that builds it; a deferral says what would bring it back.

## Next

The order of work is the first slice: tree core, the language's first slice with the query language's first forms, end
to end to SVG and a thin Canvas2D viewer, one source by region, then measure.

- **The tree core in `@molstar/core/repr`**, minimal: tree (typed columns, hashes, addresses, issues), schema, builder,
  the codec interface, JSON and S-expression codecs.
- **The language, first slice**: `spec`, `data` (a sequence, `ranges` and `values` annotations, an `offset` and a
  `codon` map), `views` (sections, tracks, `letters`, `ranges`, `values` layers); the query language's first forms.
- **End to end**: the data store, the resolved model, the resolved view, draw lists, the SVG and recording backends, a
  thin Canvas2D viewer with projected hover; views: p53 tracks, p53's protein on its gene.
- **One source by region**: FASTA with its index, or a synthetic asynchronous source first.
- **Measure**: frame time, allocation, test time, bundle size.
- **Then**: patches and `diff`, MessagePack (from io), structure maps and the Mol\* showcase, addressing and links,
  bigWig, the remaining layers and layouts, the typed configuration, the other datasets, the UX foundation stories.

## Deferred

- **The Genomes milestone**: layers across sections (ribbons, a dot plot with its own y section), synteny, contact maps,
  MAF over a genome, values for genomics (points, split, log scale, reference lines), variants by sample, matrices,
  keyed joins, gene hierarchies, a name index, two assemblies on one axis, and the formats they need (2bit, BigBed, VCF,
  PAF, GFF, chain files). The coordinate model must not rule them out.
- **Smaller layer features** (numeric filters, lane order, px heights): each with the example that needs it.
- **Definitions and uses** (`def`, arguments, templates per member of a set): when an example needs them.
- **Layouts and chrome parts as public extension points**: when an outside application needs one.
- **Settings in layers and a generated settings form**: when the typed configuration has users.
- **The compute port (workers)**: when a profile shows the main thread is the limit.
- **A JSON Schema** derived from the language description: when an outside reader or an agent tool needs to validate
  without the core.
- **Allocation guards**: for the draw and hover loops, when they exist.
