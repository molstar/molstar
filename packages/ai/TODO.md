# AI: to do

What is missing, and what was considered and deferred. One entry each, with a pointer; an entry leaves when it is done
(the spec says so then). Gathered 2026-10-16 from [\_spec](_spec/README.md) and the code. Prism's (the viewer's agent,
the assistant): [../prism/TODO.md](../prism/TODO.md).

## Decisions open

- **Reconsider zod** (David; deferred 2026-10-16). Where it is: `tools/` only among our folders (`src/tools/input.ts`: a
  field shape becomes a zod object, then the JSON Schema the model sees and the validator with its error messages;
  `src/tools/toolkit.ts`: the disclosure tools, through `zodBacked` over an explicit field shape; no zod in the public
  API). `providers/` depends on it as the AI SDK's peer, and the SDK's providers import zod themselves. Cost: ~192 kB
  minified in the assistant's chunk (English only, `scripts/esbuild/app.mjs`; ~453 kB and 95 kB gzipped with every
  locale); never in a page's startup bundle. Alternatives: JSON Schema straight from the field shapes, with a small
  validator over their closed vocabulary (string, enum, number with bounds, boolean, array, object, tagged union) and an
  open enum (listed, any string passes) for the disclosure tools, first checked by a one-off test that every viewer
  tool's schema is unchanged; JSON Schema from Mol\*'s `ParamDefinition`, which the viewer's operations already have
  (`paramsShape` maps it to field shapes today); a smaller validator. The gain is mostly code (no `toJSONSchema`
  configuration, no zod issue mapping): zod leaves the bundle only if the providers do not use the AI SDK
  ([\_spec/tools.md](_spec/tools.md) §2, "Decisions"). The question is the workspace's (one schema vocabulary beside
  `ParamDefinition` and MolViewSpec's io-ts): [../prism/TODO.md](../prism/TODO.md).

## Missing

- A failed call skips every later call in its response; the harness could skip only those that depend on it
  (`src/session/drive.ts`, the failure decision).
- More than one active turn per conversation ([\_spec/harness.md](_spec/harness.md) §5, §17).
- Persistence: a conversation that outlives the page (servers, headless hosts, a reload). The storage port went with the
  durable harness, so it is a redesign, not a backend to add; it can start from git (`e91cf9bc0`:
  `packages/ai/{harness,components,conformance}`, the host protocol, the durable core, the dropped cases)
  ([\_spec/harness.md](_spec/harness.md) §14, §17). Threads keep their row-store port: a store of the application's
  ([\_spec/threads.md](_spec/threads.md)).
- Evals: Stories' scenarios against real models, metering, baselines, the trace of each run, as a folder here (the
  viewer's first evals are its own, `prism/viewer/src/agent/eval.ts`).
- `toolscript` and `mvs-code` (code against the viewer, MVS scripting), and an MCP server over the same tools
  ([\_spec/harness.md](_spec/harness.md) §18).
- Some tests wait a fixed number of event-loop rounds (`settle()` of `@molstar/ai/testing`); waiting on what each case
  needs is the sturdier form ([\_spec/testing.md](_spec/testing.md), "Pitfalls").
- The simplification's review gate (2026-10-16; signed off by the owner the same day): CI on the pushed branch, and a
  kept record surviving a reload checked in the browser.

## Considered and deferred

- Durable nested runs: a retried parent call runs its nested runs again ([\_spec/harness.md](_spec/harness.md) §9,
  "Decisions").
