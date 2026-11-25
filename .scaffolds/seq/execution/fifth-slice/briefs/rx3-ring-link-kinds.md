# Lane RX3: kinds of links on the ring

**Base**: `packages-seq` at the orchestrator's head after RX2 merges (`git reset --hard packages-seq` first). Branch
`seq/rx3`. Model: strong. **The language, alone** if a param is needed (a minor bump; say so first in the report).

**Why**: the author, 2026-11-23, asking what a link between chromosomes means, then "so how do we differentiate between
these types of contacts in the view?", and on the orchestrator's answer (layers per kind, a colour per layer, the card
naming the kind, a story): "make it RX3 after RX2". A link's geometry says nothing of its meaning (3D contacts,
structural variants and translocations, fusion genes, synteny); the kind lives in the data, so the view shows it through
layers and fields. Read [ring-extension](../../../design/ring-extension.md) (Colour, Defaults, Configuration) and the
language's layer `filter` and `color` channels in the core `_spec`.

**Task**, a green commit at each step:

1. **A ribbons layer's colour**: `rondo:ribbons` takes a fixed colour or a scheme over a field (the layer's `color`
   channel, as other layer types read it) instead of always the block hue; the default stays the block hue
   (`rondo:hues`). Stems and stubs take the layer's colour. Kinds are layers, never mixed inside one ribbon: a bin-level
   ribbon sums its links, so colour by a field applies at the item level, and at the bins level a field colour is the
   layer's (a test of each, and of the rule's issue when a field colour meets bins).
2. **Layers per kind**: two `rondo:ribbons` layers over one source, each with its data filter on `type`, ordered (the
   rarer kind on top), widths by their own rule (count, or a field summed as `q`), each its legend entry and
   `layer.hidden` toggle; neck additivity per layer. `rondoRing` takes kinds as an option.
3. **The card names the kind**: a link's hover card lists its fields (type, the ends' strands, `q`, gene names where the
   file has them) through the `describe` hook; a bin's card its count per layer.
4. **The showcase data**: the generator's hub pairs become a few named translocations at invented coordinates (one a
   chr9–chr22 junction in the style of BCR-ABL1, named as made up), `type` `contact` and `translocation`, strands, and
   gene-name fields on the fusions; byte-identical, the BEDPE fixture regenerated, `sources.json` kept.
5. **Stories**: "Kinds of links": faint contacts in the chromosome hues, translocations in one accent over them, the
   legend toggling each; at the fit, on one chromosome, and on a pair; light and dark, canvas and SVG, headless twins.

**Where**: a worktree on `seq/rx3`. Headers (`@author David Sehnal <david.sehnal@gmail.com>`, 2026); hot files never
staged (`package.json`s, `tsconfig*.json`, lockfile, `biome.json`, `.scaffolds/seq/`, `.v6/changelogs/`, `storybook/packages.ts`,
`packages/seq/TODO.md`, `packages/seq/_measure/sizes.mjs`): list them. Anything in `@molstar/core` or `@molstar/io`
stops the lane. Sizes soft (report growth with its cause); timing and allocation bounds are gates (a frame with two
ribbons layers 0 B warm; 12,000 links in two layers under the ring's frame bound).

**Checks**: `pnpm test:unit packages/seq/ cli/seqx-render/`, `npx tsc -b packages/seq/core packages/seq/viewer
packages/seq/examples packages/seq/ext-rondo`, each `typecheck:test`, the storybook typecheck, `pnpm lint`,
`pnpm format:check`, both perf tiers to completion, `pnpm -s -C packages/seq/core size` (read the rows), `pnpm
test:medium packages/seq/`, the `.scaffolds/seq/_checks`.

**Report**: under 350 words: commits; checks with numbers; story ids; sizes; the `.v6/changelogs/` text; designs to
amend (ring-extension's Colour and Configuration).
