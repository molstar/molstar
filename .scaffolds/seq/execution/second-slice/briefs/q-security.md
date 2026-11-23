# Lane Q: the security bounds

**Goal** (plan.md, Q): `.seq/design/security.md`'s "Bounds on everything a spec sizes" table, enforced:
- every row refused with a coded issue **before allocating**, each a config item with the table's default (EC-58, 76,
  86, 91, 179), in `seqx-core` (the language's checks, the data store, the resolve passes, the sources) and the viewer's
  composition;
- a test per row at the bound and one past it (the issue, nothing allocated past the bound: check with the opt-in
  allocation tier where it matters);
- SVG output escaping (`packages/seq/core/src/render/svg.ts`) tested with hostile labels and ids (quotes, `<`, `&`,
  `]]>`, control characters, lone surrogates), and attributes;
- the headless policy's root and hosts (security.md, Defaults by host).

Not in scope: selection and emphasis (lane P edits `render` for emphasis at the same time: touch only `svg.ts`'s
escaping there), anything in `packages/core` (the asset manager's options are merged; use them).

**Read first**: `.seq/design/security.md`, `sources.md`, `data-store.md`; the core `_spec`'s Configuration table;
`packages/seq/core/src/{lang,data,resolve,context,render/svg.ts}`; `packages/seq/AGENTS.md`.

**Model**: fast (sonnet).

**Base**: `packages-seq` at `3e6f9267d`. Your worktree may be cut from master: with a clean tree, first
`git checkout -B seq/q-security 3e6f9267d`, and check `git merge-base --is-ancestor 3e6f9267d HEAD`.

**Touches**: `packages/seq/core/src/**` (not `render` beyond `svg.ts`), `packages/seq/viewer/src/{config.ts,viewer-spec.ts,mount.ts}`
where a bound is a composition setting, the `_spec`s.

**Two consumers**: one bound (the nesting depth) refused the same way by the language's check of a whole spec and by a
patch's incremental check (`checkSeqChanges`).

**Performance representative**: checking a spec at the bounds' edges stays within the whole check's time (a quick case).

**Checks**: each its own command, green only on exit 0 with a non-zero test count.
- `pnpm test:unit packages/seq/`
- `npx tsc -b packages/seq/core packages/seq/viewer packages/seq/examples`
- `pnpm -C packages/seq/core typecheck:test`, the same for `viewer`
- `pnpm lint`, `pnpm format:check`

**Headers**: a new file has the Mol\* header `Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for
more info.` with `@author David Sehnal <david.sehnal@gmail.com>`; a modified file keeps its header, its year range
extended to 2026 and that author line added if missing. No mention of `packages/sequences` or `.seq/` in code, comments,
tests or `_spec`s; code cites its own `_spec` sections.

**Fixes land with a test shown failing on the fix's parent**: the report names the command.

**Hot files, never staged**:
- `pnpm-lock.yaml`, `scripts/workspace/inventory.json`, the root `tsconfig.json`, `biome.json`, `vitest.config.ts`;
- `storybook/packages.ts`, `.v6/changelogs/` (write the changelog entry's text in your report);
- every `package.json` and `tsconfig*.json` (ask in the report: exports, `_test` nulls, dependencies);
- anything in `.seq/`.

**Stop when** the checks are green, or a check has failed twice after a fix, or the work needs a public API past the
agreed surface, or a target cannot be met without changing a design (report the measurement and the options), or a
permission is refused (report the action and why; do not work around it).

**Report**: the commits; each check with its exit code and test count; the hot-file changes you need; the changelog
entry's text; config items with their defaults (and which you were unsure of); drift (what, where, options, your
recommendation); followups.
