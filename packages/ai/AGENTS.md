# AI

The agent harness and what runs on it, the one package `@molstar/ai`. The repository's rules
([../../AGENTS.md](../../AGENTS.md)) apply; these are the package's own. Architecture and decisions:
[.v6/ai](../../.v6/ai/README.md) (`harness.md`, `tools-and-session.md`).

- **The spec** is `_spec/` (start with [_spec/README.md](_spec/README.md): the folders, how a turn binds them, where to
  start for a task). Read it before the code; a change that moves what it says updates it in the same commit.
- **What is missing** and what was deferred: [TODO.md](TODO.md) (an entry leaves when it is done).
- `src/` has the folders `messages`, `harness` (its cases are `harness/_test/`), `providers`, `tools`, `session`,
  `trace`, `threads`, `testing` (the `createSession` fixture and test doubles) and `narration`. The entries are the
  explicit `exports` of `package.json` (`_spec/README.md`, "The entries"); the root entry names what a host needs and
  nothing that pulls the AI SDK. The harness is in memory only: nothing here stores anything.
- Generic: nothing imports Prism or rxjs, nor Mol\* but `@molstar/core`'s hash file (in `messages/`); no `node:*`; zod
  only in `tools/` (Biome also allows `providers/`), and nowhere outside `packages/ai`; the AI SDK only in `providers/`.
  Biome enforces these, and which folder may import which, per folder (`biome.json`, over relative imports): a failing
  import means "not here".
- A change to a harness contract comes with a case in `src/harness/_test/`.
- Tests use the scripted model from `@molstar/ai/testing`. Anything that calls a real model costs money: ask before
  running it.

## Commands

`pnpm test:unit packages/ai` runs this package's tests with the rest of the repository's (one Vitest run);
`pnpm test:unit packages/ai/src/<folder>` one folder.
