# Lane DX: the Seq\* (next) storybook's Docs

**Base**: `packages-seq` at the orchestrator's head (`git reset --hard packages-seq` first). Branch `seq/dx`. Model:
strong. Taken early from the sixth slice (the author, 2026-11-23: "i think we can write the slice 6 docs now ?"). Beside
you: **RX3** (kinds of links: `packages/seq/ext-rondo`, `storybook/src/_stories/features/genome-ring/`,
`examples/src/genomes/ring/`). Stay out of its files; link to the ring's stories, don't edit them.

**Why** (the sixth slice's plan, "Docs"; the author, 2026-11-23): "we should write "Docs" folder in storybook in slice 6
too, that explains this very well right of the bat, so people dont expect opiniated viewer like jalview out of this. and
in general the docs should go high level around how this works etc. similar to what the docs for the previous version
do". And the principle the docs lead with ([design README](../../../design/README.md), "The engine, not the game"): "the
base library is like a "game engine", not the "game itself"".

**Read first**: the previous version's docs, for their shape, voice and live examples
(`packages/sequences/storybook/src/_stories/docs/`: `introduction.mdx`, `concepts.mdx`, `principles.mdx`, `format.mdx`,
`features.mdx`, `extending.mdx`, `query-language.mdx`, `streaming.mdx`, `color-schemes.mdx`, and `live.tsx`); the new
Seq\*'s `packages/seq/core/_spec/README.md` and `packages/seq/viewer/_spec/README.md` (the truth the docs must agree with);
`.seq/design/` (`vision.md`, `building-blocks.md`, `composition.md`, `version-1.md`, `extensibility.md`, `README.md`) for
the reasons; the existing `Seq* (next)/Design/*` pages in `packages/seq/storybook/src/_stories/design/` (those are for
contributors; the Docs are for users and hosts, and never duplicate them, they link). Do not copy the old docs: the
packages, names and language are new. `.seq/` is never cited from the docs (the `.seq` rule): cite the specs and the
stories.

**Task**: a `Docs/` folder in the Seq\* (next) storybook (`packages/seq/storybook/src/_stories/docs/`, titles
`Seq* (next)/Docs/…`, ordered), high level, each page short, with live examples (the viewer mounted on a small spec from
`seqx-examples`, as the old `live.tsx` does) and links into the stories:

1. **Introduction**: what Seq\* is and is not, first: a library that makes the pictures (the language, layers, layouts,
   navigators, data by region, hover and commands, extensions, headless and SVG), and leaves to the application how they
   are composed and navigated between; not an opinionated viewer like Jalview. One live view; the genome ring and its
   square as the worked example of an application built on the library (links to `Features/Genome Ring`).
2. **Concepts**: the building blocks: specs, data and sources, coordinate systems and maps, sections, tracks, layers,
   themes, queries, issues; how a spec becomes a picture (resolve, layout, draw lists, canvas or SVG).
3. **The language** (SeqViewSpec): its levels, the builder, versions and what a reader without an extension draws.
4. **Data**: sources and formats, read whole or by region, the store, pending and empty, the seq data port (a host's
   data through the same path).
5. **Interaction**: hover (the contract), selection, focus, navigators, commands, undo; the keyboard; tooltips and the
   card's providers.
6. **Extending**: an extension's descriptor and kinds, fallbacks, the ring as the example; card providers; what is a
   host's (the application's) choice and how it composes the parts.
7. **Embedding**: mounting the viewer in a host, the viewer spec, headless and the CLI (`seqx-render`), linking to a
   structure viewer.
8. **Performance and limits**: what it holds at which sizes, briefly, from the specs.

Each page agrees with the current `_spec`s; where a spec and a design disagree, the spec wins and the disagreement goes
in your report. Light and dark; the docs pages pass the storybook's a11y checks as the other pages do.

**Where**: a worktree on `seq/dx`. Headers on any new `.ts`/`.tsx` (`@author David Sehnal <david.sehnal@gmail.com>`,
2026); a commit per page or two; hot files never staged (`package.json`s, `tsconfig*.json`, lockfile, `biome.json`,
`.seq/`, `.v6/changelogs/`, `storybook/packages.ts`, `packages/seq/TODO.md`, `packages/seq/_measure/sizes.mjs`): if the
storybook needs registration there, list it. No code change outside the storybook's docs folder; a doc that finds a bug
reports it.

**Checks**: the storybook typecheck, `pnpm test:unit packages/seq/storybook`, `pnpm lint`, `pnpm format:check`, and a
storybook build or dev start that renders every new page without errors (say which you ran and how you checked).

**Report**: under 300 words: commits; the pages and their story ids; checks; spec–design disagreements found; hot-file
needs; the `.v6/changelogs/` text.
