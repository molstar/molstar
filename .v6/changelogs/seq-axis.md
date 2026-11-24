## Seq* (next): the axis: direction, collapsed introns, letters through a path

- A section's axis can be reversed or follow a path's direction (`direction`: `genomic`, `reversed`, `follow`,
  `{ along }`), and show only some of its positions (`pieces`, `collapsed`: `coding`, ranges or a map, with `pad`);
  breaks are drawn as hairlines, hit and lit.
- Letters read through one-to-one and codon paths: a residue a codon wide, three-letter codes, the complement on a
  reversed strand, a codon conflict marked under its residue's cell, carets and gap marks.
- Viewer commands `section.reverse` and `section.collapse` (Shift+R, Shift+C); `selection.acrossBreaks`;
  `axis.directionMark`; `renderSvg` accepts a viewport in positions and per-section viewports.
