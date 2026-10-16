# AI: to do

What is missing, and what was considered and deferred. One entry each, with a pointer; an entry leaves when it is done
(its design or plan says so then). Gathered 2026-10-16 from [.v6/ai](../../.v6/ai/README.md), [\_spec](_spec/README.md)
and the code. Prism's (the viewer's agent, the assistant): [../prism/TODO.md](../prism/TODO.md).

## Decisions open

- **Reconsider zod** (David). Where it is: `ai-tools` only among our packages (`tools/src/input.ts`: a field shape
  becomes a zod object, then the JSON Schema the model sees and the validator with its error messages;
  `tools/src/toolkit.ts`: the disclosure tools, through `zodInput`; `tools/src/zod.ts`: `zodInput` for hosts that
  declare tools in zod). `ai-model-adapters` depends on it as the AI SDK's peer, and the SDK's providers import zod
  themselves. Cost: ~192 kB minified in the assistant's chunk (English only, `scripts/esbuild/app.mjs`; ~453 kB and 95
  kB gzipped with every locale); never in a page's startup bundle. Alternatives: JSON Schema straight from the field
  shapes, with a small validator over their closed vocabulary (string, enum, number with bounds, boolean, array, object,
  tagged union), keeping `zodInput` as an opt-in subpath; JSON Schema from Mol\*'s `ParamDefinition`, which the viewer's
  operations already have (`paramsShape` maps it to field shapes today); a smaller validator. zod leaves the bundle only
  if the adapters do not use the AI SDK ([tools-and-session.md §1–2](../../.v6/ai/tools-and-session.md)).

## Missing

- A failed call skips every later call in its response; the harness could skip only those that depend on it
  (`session/src/drive.ts`, the failure decision).
- More than one active turn per conversation ([harness.md §5](../../.v6/ai/harness.md);
  [\_spec/README.md](_spec/README.md), "State of the code").
- SQLite storage for the harness and fibers, for servers and headless hosts (harness.md §14,
  [plan.md](../../.v6/ai/plan.md) "Later", [\_spec/fibers.md](_spec/fibers.md)).
- `@molstar/ai-evals`: Stories' scenarios against real models, metering, baselines, the trace of each run (plan.md
  "Later"; the viewer's first evals are its own, `prism/viewer/src/agent/eval.ts`).
- `ai-toolscript` and `mvs-code` (code against the viewer, MVS scripting), and an MCP server over the same tools
  (plan.md §2.6, [analysis/viewer-mode.md §7](../../.v6/ai/analysis/viewer-mode.md)).
- The conformance suite waits fixed event-loop rounds in places and retries a file; waiting on what each case needs
  removes the retries (plan.md "Later").

## Considered and deferred

- Durable nested runs: a nested run re-runs instead ([harness.md §9](../../.v6/ai/harness.md), README decision 5).
- Data migration: a turn written by another schema version is stopped or reset ([README](../../.v6/ai/README.md)
  decision 9).
