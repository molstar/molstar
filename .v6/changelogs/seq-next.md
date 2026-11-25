## Seq\* (next): a new implementation of SeqViewSpec and its viewer

- `@molstar/seqx-core`, `@molstar/seqx-viewer` and `@molstar/seqx-examples` (`packages/seq`, temporary names until
  they replace the current Seq\* packages): SeqViewSpec, a language for sequence views as a tree (JSON, printed as
  S-expressions), over one coordinate model where data on one system shows on every system it maps to.
- The core: the language and its builder; queries (ranges, items, nodes, labels in `index` and `ordinal` numberings);
  a source protocol read region by region (tiles, levels, summaries) with a byte-bounded data store; FASTA with its
  index; dynamic nodes and track templates; resolve through codon and offset maps; draw lists for letters, ranges and
  values; SVG export (`renderSvg`), headless.
- The viewer: `createSeqViewer` composed from an explicit spec (`DefaultSeqViewerSpec`), DOM-free models (spec, view
  with scroll groups and follow, hover projected across maps, store), a Canvas2D backend on three stacked canvases,
  themes from the `@molstar/ui` root, an asset policy over the asset manager.
- Performance, headless on the 3-million-position fixtures: pan and zoom frames under 2.3 ms; drawing, hover, selection
  and cursor redraws allocate nothing once warm; bundles 133 kB (reading a spec), 229 kB (headless SVG), 318 kB (the
  viewer), minified.
- Storybook: `Seq* (next)`, with the p53 examples (Tracks, Protein on its gene, Gene by region).
