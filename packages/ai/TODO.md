# AI: to do

What is missing, and what was considered and deferred. One entry each, with a pointer; an entry leaves when it is done
(its design or plan says so then). Gathered 2026-10-16 from [.v6/ai](../../.v6/ai/README.md), [\_spec](_spec/README.md)
and the code. Prism's (the viewer's agent, the assistant): [../prism/TODO.md](../prism/TODO.md).

## Decisions open

- **Reconsider zod** (David; slice 12b of the simplification, deferred:
  [tool-schemas.md §3](../../.v6/ai/simplification/tool-schemas.md)). Where it is: `tools/` only among our folders
  (`src/tools/input.ts`: a field shape becomes a zod object, then the JSON Schema the model sees and the validator with
  its error messages; `src/tools/toolkit.ts`: the disclosure tools, through `zodBacked` over an explicit field shape; no
  zod in the public API). `providers/` depends on it as the AI SDK's peer, and the SDK's providers import zod
  themselves. Cost: ~192 kB minified in the assistant's chunk (English only, `scripts/esbuild/app.mjs`; ~453 kB and 95
  kB gzipped with every locale); never in a page's startup bundle. Alternatives: JSON Schema straight from the field
  shapes, with a small validator over their closed vocabulary (string, enum, number with bounds, boolean, array, object,
  tagged union); JSON Schema from Mol\*'s `ParamDefinition`, which the viewer's operations already have (`paramsShape`
  maps it to field shapes today); a smaller validator. zod leaves the bundle only if the providers do not use the AI SDK
  ([tools-and-session.md §1–2](../../.v6/ai/tools-and-session.md)). The question is the workspace's (one schema
  vocabulary beside `ParamDefinition` and MolViewSpec's io-ts): [../prism/TODO.md](../prism/TODO.md).

## Missing

- A failed call skips every later call in its response; the harness could skip only those that depend on it
  (`src/session/drive.ts`, the failure decision).
- More than one active turn per conversation ([harness.md §5](../../.v6/ai/harness.md);
  [\_spec/README.md](_spec/README.md), "State of the code").
- Persistence: a conversation that outlives the page (servers, headless hosts, a reload). The storage port went with the
  durable harness, so it is a redesign, not a backend to add ([harness.md §17](../../.v6/ai/harness.md),
  [plan.md](../../.v6/ai/plan.md) "Later"). Threads keep their row-store port: a store of the application's
  ([\_spec/threads.md](_spec/threads.md)).
- Evals: Stories' scenarios against real models, metering, baselines, the trace of each run (plan.md "Later"; the
  viewer's first evals are its own, `prism/viewer/src/agent/eval.ts`).
- `toolscript` and `mvs-code` (code against the viewer, MVS scripting), and an MCP server over the same tools (plan.md
  §2.6, [analysis/viewer-mode.md §7](../../.v6/ai/analysis/viewer-mode.md)).
- Some tests wait a fixed number of event-loop rounds (`settle()` of `@molstar/ai/testing`); waiting on what each case
  needs is the sturdier form ([\_spec/testing.md](_spec/testing.md), "Pitfalls").

## Considered and deferred

- Durable nested runs: a retried parent call runs its nested runs again ([harness.md §9](../../.v6/ai/harness.md),
  [README](../../.v6/ai/README.md) decision 5).
