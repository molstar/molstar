# AI: to do

The decisions still open. One entry each, with a pointer; an entry leaves when it is decided (the spec says so then).
Gathered 2026-10-16 from [\_spec](_spec/README.md) and the code. Prism's (the viewer's agent, the assistant):
[../prism/TODO.md](../prism/TODO.md).

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
