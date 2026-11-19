# Lane GR: research, a genome ring like Rondo's (read-only)

**Question** (the author, 2026-11-19, with three screenshots of https://www.rondo.app/): "for the genome in step 5,
would we be able to have a view like this or close to it? … maybe run a research for this … we can do something like
this with our representations, perhaps even as an extension to show this is "easy" with the current system?"

**What the screenshots show** (describe them back in your report as you understand them):
1. **The genome ring**: chromosomes 1–22 and X as blocks on a ring with gaps and a gap at the top, each labelled
   outside, each with a cytoband ideogram strip at its outer edge, colored per chromosome (a hue wheel). Crosslinks
   (pairs of positions between or within chromosomes, about 12,000 from one chromosome alone) drawn inside as curves,
   **bundled**: from each chromosome block a filled "trumpet" fans into its links, links between two chromosomes merge
   into thick ribbons, a link within a chromosome a loop back to it.
2. **Hover on a chromosome** (3): its links highlighted, the rest faded; Mb ticks around it; a tooltip
   "crosslinks: 11972, 0-198022430"; the blocks show finer fans per band.
3. **One chromosome as its own ring** (chromosome 4): its cytobands as the ring, Mb ticks outside, the links within it
   as ribbons inside the ring between bundled "bulbs" per region, the links to other chromosomes leaving outward as
   thin colored lines.
4. **A navigator**: a rounded square of the chromosomes as colored cells with their numbers (nested, 1 and X at the
   centre), the current chromosome marked.

**Answer, in order**:
1. **Rondo itself**: what it is, its data (crosslinks of what: Hi-C, XL-MS?), how it bundles (hierarchical edge
   bundling? aggregation by bins?), its levels of detail and interactions. From its public pages only; no sign-in, no
   download. If the site gives little, say so and work from the screenshots.
2. **Against Seq\*'s designs and code** (`.seq/design/`: layouts.md (ring, blocks), grids-and-rows.md and
   `pairs`/arcs, navigators.md, themes.md, extensibility.md, sources.md and data-store.md (by region, levels),
   performance.md, features.md's deferred rows: layers across sections EC-209, two assemblies EC-210, contact maps
   EC-223; and `packages/seq/core/src`): for each element above, **built** (the module and a story that shows it,
   e.g. Complex's contacts as arcs across blocks, the Plasmid ring), **close** (what small addition), or **missing**
   (what design it needs). Particularly: arcs between blocks on a ring at genome scale (12,000+ pairs: culled? summed
   by level?); bundling (a pairs reduction by bin pairs at the level that fits, drawn as ribbons whose width is the
   count: is that a new layer type, a theme, or a `pairs` level?); the trumpets and bulbs (a shape per bin's fan); a
   section that is one block of another (zoom from the genome ring into chromosome 4's own ring, its outward links);
   the hue-wheel coloring per block; the square navigator.
3. **As an extension**: could a package outside `packages/seq` (an `extensions/`-style package, registry entries only:
   a layer type, a reduction, a navigator part) build it on the public API, by extensibility.md's rules? What would it
   need from core that is not public today (each a finding)? A sketch: the spec a host writes, the registry entries,
   the data (a generated crosslink set over GRCh38's lengths and cytobands; cytobands from UCSC are public, licence to
   check), roughly how much code.
4. **Recommendation**: what is the smallest demo that shows "easy" honestly; which slice (the fifth, Genomes, by the
   author's decision; or an extension demo earlier); what it would teach the designs.

**Model**: fast (sonnet). **Where**: the main checkout, read only: no file written, no installs or builds, no git
writes. Web: Rondo's public pages and papers or docs about it, nothing that needs an account.

**Report**: the four answers, each claim about Seq\* with its file and line; a table of the elements (built, close,
missing); open questions for the author.
