# State: adopt-principles (2026-11-17)

## Done

Every step of plan.md. The instance in `.seq/`, the package in `.seq/oz/`, the checks wired and shown to fail on
planted errors, the sources and root designs, `packages/seq/{core,viewer,storybook}` registered (pnpm workspace,
inventory, root tsconfig references, Storybook group `Seq* (next)`), `.seq/` excluded from Biome and Prettier,
`.seq/tmp/` gitignored.

## Verify

- The three checks (`.seq/AGENTS.md`, Checks): ok.
- `npx tsc -b packages/seq/core packages/seq/viewer`, `npx tsc -p packages/seq/storybook/tsconfig.json`: ok.
- `pnpm test:unit packages/seq`: passes (the filter also matches `packages/sequences`).
- `pnpm lint`, `pnpm format:check`: ok.
- `pnpm check:workspace`: fails on this checkout before and after the change, on exports of unbuilt packages
  (`lib/` missing); no message names the new packages. Needs `pnpm build` first to read clean.
- Storybook: the `Seq* (next)` group shows Introduction, Building blocks, Format, Sources, Examples.

## Decisions taken

In `.seq/sources/author/README.md` and the designs. One not asked: the package lives in `.seq/oz/`, not in
`.seq/principles/` as meso\* has it, because the package's checks require the package and the instance apart.

## Next

The language description (`.seq/status.md`, Open 2), design first.

Closed 2026-11-17, recorded in `.seq/records/adopt-principles/2026-11-17-scaffold.md`.
