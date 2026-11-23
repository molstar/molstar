# Plan: adopt the principles in .seq and scaffold packages/seq

The adopting procedure (`.scaffolds/principles/_core/patterns/adopting.md`) plus the author's request for the initial scaffold.

## Checklist

- [x] 0. Ask first: the author's answers in `.scaffolds/seq/sources/author/README.md` (format, sources, the rest as recommended).
- [x] 1. Instance: manifest at `.scaffolds/seq/principles/README.md`; the package copied to `.scaffolds/principles/`; meso\*'s queries and
      agent-interface as the instance's own modules; provenance and adaptations recorded.
- [x] 2. Coexist: `.scaffolds/seq/AGENTS.md`; the one pointer from `packages/seq/AGENTS.md` (no pointer from the root
      `AGENTS.md`: the scaffold is temporary).
- [x] 3. Skills: the instance skill at `.scaffolds/seq/_skills/seq/`; skills read as files, not linked (manifest, Folder map).
- [x] 4. Checks: id, path, merge, with `--allow` limiting them to the scaffold; each fails on a planted error.
- [x] 5. Sources and the root designs: vision, building blocks (provisional), format, sources, packages, examples.
- [x] 6. `packages/seq/{core,viewer,storybook}`: package.json, tsconfig, inventory, root references, `_spec` stubs,
      `TODO.md`, the Storybook entry `Seq* (next)` with the docs pages, one smoke test per package.
- [x] 7. Gate (`check:workspace` reads clean only after `pnpm build`: state.md, Verify): `pnpm check:workspace`, `pnpm lint`, `pnpm format:check`, `npx tsc -b` for the new projects, the quick
      lane for `packages/seq`.
- [x] 8. Record: `.scaffolds/seq/records/adopt-principles/`.

## Testing strategy

The scaffold has no language code: a test per package that its entry imports and its folder graph holds; the Storybook
docs pages build in the Storybook typecheck.
