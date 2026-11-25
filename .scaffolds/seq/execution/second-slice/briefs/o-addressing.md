# Lane O: addressing in Seq\*

**Goal**: build `.scaffolds/seq/design/addressing.md` in `@molstar/seqx-core` and the viewer's commands (plan.md, O):
- **ids and addresses to what they name and back**: an address (core's `parseAddress`/`selectAddress` grammar, as is)
  resolving to a spec node, a node data makes (forms), or a Location (ranges, rects) through the query model; and
  from a Location or node to its address, the shortest stable form (a ref if one names it, else ids and forms);
  bounded (limit, depth), more matches an issue, never truncated; an address that no longer resolves reported with
  what it named;
- **user-defined refs and bookmarks**: the `refs` container and the `ref` kind in `SeqLanguage` (`to`, `color`,
  `label`; a set of addresses for a location set), `/refs/<id>` usable wherever an address is; adding a bookmark is
  one patch (`applyPatch` from `@molstar/core/repr/patch`); a ref to a dynamic node waits (an issue, not an error)
  while its node is not held, and resolves when it arrives;
- **actions and links**: commands on addresses (focus, reveal, select, highlight on a named channel, hover) in the
  viewer, so a link, a key, a test and an agent drive the view the same way; links in Mol\*'s link grammar
  (`!seq-focus=/refs/dbd`, percent-encoded as a whole), parsed and printed both ways, no markup from data; copy link;
- **addresses of the session state**: the selection and the focus readable and settable as addresses.

Not in scope: editing through `SpecModel` and history (lane N); selection's behavior and highlights' drawing beyond
what a command needs to reach them (lane P); markdown rendering.

**Agreed surface**: nothing new in `@molstar/core`. If core's address grammar or `selectAddress` lacks what the design
needs, stop and report it.

**Read first**: `.scaffolds/seq/design/addressing.md`, `patches.md` (refs following moves, dangling and rebinding), `changes.md`,
`viewer-models.md`, `composition.md`, `language-description.md`; `packages/seq/core/_spec/README.md` (Language, Query),
the viewer `_spec`; `packages/core/src/repr/tree/address.ts`; `packages/seq/core/src/{lang,query,resolve}` and
`packages/seq/viewer/src/{commands.ts,models}`; `packages/seq/AGENTS.md`.

**Model**: strong (opus).

**Base**: `packages-seq` at `595c6a26d`. Branch `seq/o-addressing`; commit on it only, at each small green step.

**Touches**: `packages/seq/core/src/**` (a new `address/` folder or as the core `_spec` places it; `lang` for `refs`),
`packages/seq/viewer/src/**` (commands, models, ports), `packages/seq/examples/src/**`, `packages/seq/storybook/**`
stories under `Seq* (next)/Features/Addressing` and their headless twins; their `_spec`s. Not `packages/core`.

**What is drawn**: a bookmarks story (refs listed, each focusing the view; adding one is one patch) and an addressing
story (a short text with links that focus, highlight and select; one location in two numberings; copy link), each with
its headless twin composed as the story composes it, on a fake clock.

**Two consumers**: the viewer and an export (SVG) resolving the same address in one test, to the same Location.

**Performance representative**: resolving 1,000 addresses (ids, forms and refs mixed) on the p53 spec: a quick case, and
a heavy bound in `itWithin` (propose it from measurement).

**Configurable over fixed**: the selection bounds' defaults, the link key prefix (`seq-`), the copy-link preference
order, in the `_spec`s' Configuration tables.

**Checks**: each its own command, green only on exit 0 with a non-zero test count.
- `pnpm test:unit packages/seq/`
- `pnpm test:heavy packages/seq/`
- `npx tsc -b packages/seq/core packages/seq/viewer packages/seq/examples`
- `pnpm -C packages/seq/core typecheck:test`, the same for `viewer` and `examples`, and `pnpm -C packages/seq/storybook typecheck`
- `pnpm lint`, `pnpm format:check`

**Headers**: a new file has the Mol\* header `Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for
more info.` with `@author David Sehnal <david.sehnal@gmail.com>`; a modified file keeps its header, its year range
extended to 2026 and that author line added if missing. No mention of `packages/sequences` or `.scaffolds/seq/` in code, comments,
tests or `_spec`s; code cites its own `_spec` sections.

**Fixes land with a test shown failing on the fix's parent**: the report names the command.

**Hot files, never staged**:
- `pnpm-lock.yaml`, `scripts/workspace/inventory.json`, the root `tsconfig.json`, `biome.json`, `vitest.config.ts`;
- `storybook/packages.ts`, `.v6/changelogs/` (write the changelog entry's text in your report);
- every `package.json` and `tsconfig*.json` (ask in the report: exports, `_test` nulls, dependencies);
- anything in `.scaffolds/seq/`.

**Stop when** the checks are green, or a check has failed twice after a fix, or the work needs a public API past the
agreed surface, or a target cannot be met without changing a design (report the measurement and the options), or a
permission is refused (report the action and why; do not work around it).

**Report**: the commits; each check with its exit code and test count; the hot-file changes you need; the changelog
entry's text; config items with their defaults (and which you were unsure of); drift (what, where, options, your
recommendation); followups.
