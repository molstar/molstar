# Plan

## Where it comes from

The harness was rewritten in Stories Untold (`lib/anystar`, branch `harness-next`, built 2026-10-05) against a frozen
spec, conformance first, and Stories ran on it from then on. On 2026-10-11 the five packages moved here unchanged in
logic, as `@molstar/ai-*`; Stories' own code (its tools, drive and agent) stayed behind. What carried over from that
work:

- **Conformance first.** A contract change comes with a conformance case, never a quiet edit; the suite runs on every
  backend.
- **One home per concern** ([harness.md](harness.md) §1). A change that moves a concern to another home is the thing to
  look for in review.
- **Import rules fail the build** (`biome.json`): no package imports Mol\*, Prism or rxjs; `ai-agent` imports nothing of
  the workspace; `ai-harness` only `ai-agent`; no `node:*` outside `ai-components`. The original stack had drifted one
  reasonable local change at a time until concerns had no clear home; a rule that fails says "not here".

## 1. The move (2026-10-11)

Done:

- **The workspace**: `packages/ai/*`, the catalog (the AI SDK pinned to its major, zod as its peer), Node 22.13 or later
  everywhere (`node:sqlite`), the import rules (`biome.json`), prism-ui's opt-in boundary.
- **The code**, logic unchanged: relative imports to `.js`, the header, the repository's formatting and lint, `lib`
  ES2024, Vitest to Jest (the suite's `expect(value, message)` became `expectWith`), the old spec's citations pointed at
  [harness.md](harness.md)'s sections.
- **The tests run apart** (`pnpm jest:ai`, `packages/ai/jest.config.cjs`): the AI SDK is ESM only and the repository's
  Jest is CommonJS, so that config transforms it. It goes when Jest moves to ESM (`.v6/plans/checklist.md`).
- **The docs**: this folder, short; the spec, its decisions log and the rewrite's plan stay in Stories Untold's archive.

## 2. Next

In this order (reordered 2026-10-12; [analysis/viewer-mode.md](analysis/viewer-mode.md) §7):

1. **The tool layer and the session**: `@molstar/ai-tools` (done) and `@molstar/ai-session`, ported from Stories' agent
   with its story cut out ([tools-and-session.md](tools-and-session.md)).
2. **`@molstar/ai-trace` and the Inspector** (Trace, Tools, Limits), on mock data until a session feeds them.
3. **`@molstar/ai-fibers`**: the kernel, the host, the backends in `ai-components`, the conformance.
4. **PD → field shape** (and so JSON Schema), for commands as tools.
5. **Pause**: one plan with the Prism UI's, and the sequence viewer. Then Prism hosts it (the integration step of
   [../prism-ui/plan.md](../prism-ui/plan.md): commands as tools, `execute` as the handler, the drive, origins, one undo
   step per turn through the mode's checkpointer, [harness.md](harness.md) §18), the skills, `ai-toolscript` and
   `mvs-code`, the MCP server over the same tools.

Later:

- **Evals**: `@molstar/ai-evals`, the common part of Stories' evals (`evals/`, `apps/eval`: scenarios run against real
  models, metering, baselines, the trace of each run), with each mode's scenarios in its `agent` folder.
- **Jest to ESM** for the whole repository, then `jest:ai` folds into `jest`.
- **The suite waits on conditions.** Some cases wait a fixed number of event-loop rounds for work that includes
  WebCrypto hashing; on fake IndexedDB under load that loses the race now and then (it did under Vitest too), so that
  file retries its cases. Waiting on what each case needs (as S3 now waits for its tool to start) removes the retries.
