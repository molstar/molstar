# AI

The agent harness and what runs on it. The repository's rules ([../../AGENTS.md](../../AGENTS.md)) apply; these are the
AI packages' own. Architecture and decisions: [.v6/ai](../../.v6/ai/README.md) (`harness.md`, `tools-and-session.md`).

- `agent`, `harness`, `components` (platform code: storage, the supervisor), `model-adapters`, `conformance` (the suite
  every backend runs, and test doubles), `tools`, `session`, `trace`, `fibers`.
- Generic: no package imports Mol\*, Prism or rxjs; `node:*` only in `ai-components`; zod only in `ai-tools`. Biome
  enforces these per package.
- A change to a harness contract comes with a conformance case.
- Tests use the scripted model from `@molstar/ai-conformance/testing`. Anything that calls a real model costs money: ask
  before running it.

## Commands

`pnpm test:unit packages/ai` runs these packages' tests with the rest of the repository's (one Vitest run);
`pnpm test:unit packages/ai/<pkg>` one package.
