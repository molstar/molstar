## Seq* (next): tooltips, and the hover card composed by providers

- A `tooltip` node adds rows to the hover card on drawn marks: by data (an item's own fields) or by range from a table (inline rows or a source's ranges, read by region), as f-string templates with inline Markdown and in-view `!seq-` links (language 0.9.0).
- The hover card is composed by card providers registered in the viewer spec (the viewer's own rows are built-in providers; a host adds or replaces one).
- A view's `export.titles` titles SVG items with their names or their cards.
- The card caps its width and wraps long values, can be entered when it holds a link, and titles a long VCF REF by its length.
