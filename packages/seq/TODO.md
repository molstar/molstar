# TODO

What is not built. An item leaves in the commit that builds it; a deferral says what would bring it back.

## Next

- **The tree core in `@molstar/core`** (tree, schema as data, patches, addresses, builder engine), which the language
  description builds on.

- **The language description**: kinds, where each may sit, ids, params with types, defaults, bounds and docs, refs;
  validation reporting every issue; defaults; the builder's types; the S-expression printer and reader. First slice:
  `spec`, `data` with `sequence` and `annotation` (ranges), one `view` with sections, tracks and the `letters` and
  `ranges` layers.
- **Patches**: `set`, `unset`, `insert`, `remove`, `move` by address; batches with inverses; `diff`.
- **The source protocol**: coverage and region requests answered as the model's tables; inline data and FASTA by region
  first.
- **The coordinate model**: systems, maps as segment tables with a direction on every segment, paths, projection.
- **The render path**: the resolved view and draw lists as public forms, the SVG and recording backends, `renderSvg` and
  a render CLI, headless in Node; the Canvas2D backend in the viewer.
- **The first example end to end**: the p53 text view, resolved, laid out, drawn to SVG and Canvas2D, in Storybook.
- **Configuration and extension points**: the settings layers with effective-value lookup and a generated settings form;
  the descriptor shape; the first point built as an extension (layer types), with a worked example in Storybook.
- **UX foundations stories**: color, type, grid, density, shape, emphasis, motion, each drawn on real data beside its
  violation; the interaction grammar's conformance tests.
- **Structure maps**: the structure system and the `structure` map kind; a SIFTS-style residue table source.
- **The Nightingale reference adapter** (`packages/seq/nightingale`): a subset of the resolved view, in Storybook and a
  test.
- **Addressing**: the address forms, parsing and printing, commands taking addresses, Mol\*-grammar links, copy link.
- **bigWig**: the first genome format runtime, by region and zoom level.
- **Showcases** (storybook package): services composed by an application; predictions onto a structure through mmCIF
  `label` and `auth` numbering with Mol\*; addressing.
- **Datasets**: p53, globins, antibody, plasmid, tRNA, two toy genomes; many views over each, editing views included.

## Deferred

- **Layers across sections** (ribbons between two sections, a dot plot with its own y section): when an example needs
  them; the coordinate model must not rule them out.
- **More genome format runtimes** (2bit, BigBed, VCF with tabix, PAF): after bigWig, each with the dataset that needs
  it.
- **Smaller layer features** (`split` on values, log scale, points, numeric filters, lane order, joins to keyed tables,
  px heights): each with the example that needs it.
- **A JSON Schema** derived from the language description: when an outside reader or an agent tool needs to validate
  without the core.
- **Allocation guards**: for the draw and hover loops, when they exist.
