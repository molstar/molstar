# Tool schemas: zod isolated in `packages/ai`

**The goal** (David, 2026-10-16): zod stays inside `packages/ai`; no zod schemas anywhere else in the app. Tools are
declared outside `packages/ai` only as field shapes (or Mol\*'s `ParamDefinition`, through `paramsShape`).

That is true of the code today (§1) but not enforced: nothing outside `packages/ai` imports zod, yet no rule stops it,
and `@molstar/ai-tools/zod` (`zodInput`) invites a host to declare tools in zod. The isolation needs three small
changes (§0); replacing zod inside `packages/ai` with our own schema code (§3) is a separate, optional simplification.

## 0. The isolation

1. **No zod in the public API**: drop the `./zod` entry (`zodInput`, `verbatim`; no importer, dead-code.md §1) and keep
   no zod type in an exported signature. `ToolInput` (`tools/src/input.ts:19–25`) already names none ("zod is this
   package's own business, so no type here names it").
2. **One lint rule for the repository**: a Biome override on every source file except `packages/ai/**` that restricts
   `zod` and `zod/**` ("zod is the AI package's own: declare tool inputs as field shapes"). Today zod is fenced only
   inside five of the AI blocks (`biome.json:135, 185, 229, 319, 415`); outside `packages/ai` the only guard is
   `pnpm check:workspace` refusing an undeclared import (`scripts/workspace/check.mjs:226`), which a package can lift by
   adding zod to its `package.json`. Biome's negated `includes` (`"!packages/ai/**"`) already appear in this file, so
   the override is expressible.
3. **Inside `packages/ai`**, zod where it is needed: `providers/` (the AI SDK's peer) and `tools/` (as long as the
   schema bridge uses it). The folder rules (packaging.md §5) forbid it in the other folders.

Size: a few lines of `biome.json`, the `./zod` entry, and `zod.ts` (118) with its test (205). Sonnet can do it (plan.md
slice 12a).

The rest of this file is the optional step: our code could drop zod entirely, since it already has the schema layer
(field shapes) and uses zod only as its back end. zod would still be installed, because the AI SDK needs it, and still
be in any bundle that calls a model.

## 1. Where zod is today

| File (`packages/ai/tools/src/`) | Lines | What zod does there                                                                                  |
| ------------------------------- | ----: | ---------------------------------------------------------------------------------------------------- |
| `input.ts`                      |   164 | a field shape becomes a zod object (`fieldSchema`, `objectSchema`); `z.toJSONSchema` with a config gives the JSON Schema the model sees (`JSON_SCHEMA`, which also strips the safe-integer bounds zod adds); `safeParse` checks the arguments; `argumentsMessage` turns the first `ZodError` issue into one message the model can act on |
| `toolkit.ts:400–404, 488–491`   |     — | the two disclosure tools (`disclose_skills`, `disclose_tools`) declare `{ skills: string[] }` / `{ tools: string[] }` in zod, with the ids as an enum **in the schema only** (`.meta(enumOf(…))`), so an unknown id still parses and gets an answer |
| `zod.ts`                        |   118 | `zodInput` and `verbatim`: a host may declare a tool's input in zod; `shapeOf` converts the zod object back to a field shape (for normalizing) by reading zod's private `_zod.def` |

No package outside `packages/ai/tools` imports zod (`git grep "from 'zod'"`), and nothing imports `@molstar/ai-tools/zod`
(dead-code.md §1). Prism declares every tool as a field shape: `packages/prism/viewer/src/agent/params.ts` maps a
command's `ParamDefinition` to field shapes (`paramsShape`). `packages/prism/viewer/src/agent/skills.ts:10` and
`packages/prism/assistant/src/spec.ts:13` keep zod and the AI SDK out of the page's startup bundle by lazy imports.

So the "simple tool schema definition layer" exists: `FieldShape` (`shape.ts`: string with enum / format /
verbatim, number and integer with bounds, boolean, array with a max, object, tagged variant, json; optional, default,
description), its types (`ArgsOf`), and `normalize` (trim, coerce, case-fold enums). zod is only how it is turned into
JSON Schema and checked.

## 2. zod and the AI SDK

The SDK packages declare zod as a peer (`ai@7.0.137`, `@ai-sdk/{anthropic,openai,google}@4`, `openai-compatible@3`,
`provider-utils@5`: `"zod": "^3.25.76 || ^4.1.8"`) and import it at module level to parse provider responses (`ai` 8
imports of `zod/v4`; `@ai-sdk/anthropic` 61; `@ai-sdk/provider-utils` 91 of `zod/v3`, 6 of `zod/v4`, 7 of
`zod/v4/core`, counted in their `dist`). The adapters already pass tool schemas to the SDK as plain JSON Schema
(`jsonSchema` from `ai`, `model-adapters/src/common/stream.ts:20`), so no zod object crosses that boundary.

Consequences:

- `zod` stays a dependency of `@molstar/ai`, as the SDK's peer, whatever the tools do.
- A page that calls a model loads zod with its first provider (the lazy `./providers/*` entries, packaging.md §2). The
  ~192 kB that `packages/ai/TODO.md` ("Reconsider zod") counts in the assistant's chunk moves to whichever chunk loads
  first, so the saving is only for a page that builds the agent and never calls a model (inferred: esbuild puts a
  module shared by lazy chunks in a shared chunk; not measured here).
- zod leaves the bundle only if the providers stop using the AI SDK, which the target keeps (providers.md).

The gain is therefore mostly code: one layer of ours instead of ours on top of zod, no reading of zod's private
`_zod.def` (`zod.ts`), no `toJSONSchema` configuration to keep in step with zod's releases, and messages built from our
own checks instead of mapped from zod's issue codes (`input.ts:107–127`).

## 3. The layer, proposed

In `src/tools/input.ts` (renamed `schema.ts` if wanted), over the existing `FieldShape`:

- **`jsonSchemaOf(shape)`**: draft 2020-12, input side, as today's output: `type`, `properties`, `required` (fields
  neither optional nor defaulted), `enum`, `minimum` / `maximum` (only the bounds given), `items` / `maxItems`, `oneOf`
  of objects with a `const` tag for a variant, `{}` for json, `default`, `description` (with "Format: …" appended).
  Inferred ~60 lines.
- **`check(shape, args)`**: walks the shape once, collecting the typed value (defaults applied, unknown keys dropped, as
  zod's non-strict object does) or the first failure as `{ path, expected, got }`; the message text is today's
  (`'mark.type: expected one of box, arrow (got "cone").'`, `'mark: required (What to draw.).'`). Inferred ~120 lines.
- **An open enum**: `{ kind: 'string', enum, open: true }` (name to choose): the values are listed in the schema but any
  string passes. That is what the disclosure tools need (`toolkit.ts:400–404`); they become field shapes.
- **`ToolInput` stays** (`{ shape, jsonSchema, parse }`, `input.ts:19–25`): a host that wants zod (or anything else)
  builds a `ToolInput` itself. That keeps the escape hatch without zod in our code.

Goes: `zod.ts` (118) and `_test/zod.test.ts` (205), the `./zod` entry, the zod imports in `input.ts` and `toolkit.ts`,
`zodBacked`, `JSON_SCHEMA`. Net: roughly 164 + 118 + 205 = 487 lines replaced by ~180 source lines (inferred).

**Checks that keep the behavior**:

- `tools/src/_test/input.test.ts` (169 lines, 8 cases) pins the JSON Schema and the model-facing messages for every
  kind; it should pass unchanged against the new layer.
- Before removing zod, a one-off test generates today's JSON Schema for every tool the viewer registers
  (`createViewerToolkit`, `packages/prism/viewer/src/agent/toolkit.ts`) and the disclosure tools, and asserts the new
  generator's output is identical, so the models see the same tool definitions (and the providers' cached prefixes do
  not change).
- `normalize.test.ts`, `tool.test.ts`, `toolkit.test.ts` and the viewer's `agent/_test/params.test.ts` run as they are.

## 4. What breaks outside

Nothing imports zod or `./zod` outside `packages/ai/tools`. The comments that name zod as a reason for lazy loading
(`packages/prism/viewer/src/agent/skills.ts:10`, `_test/skills.test.ts:16`, `packages/prism/viewer/src/spec.ts:124`,
`packages/prism/assistant/src/spec.ts:13`) stay true (the AI SDK still brings it) but can say "the AI SDK" instead.
Docs: `.v6/ai/tools-and-session.md` §1–2 (input "as field shapes (or zod)"), `packages/ai/_spec/tools.md:12,52`, the
Biome rule "zod is ai-tools' own" (it becomes "zod: providers only, as the SDK's peer"), `packages/ai/TODO.md`'s
"Reconsider zod" entry (another session's file: answered by this, removed in slice 10).

## 5. Where it goes in the plan

- **12a, the isolation** (§0): after the collapse (slice 6), or before it on today's packages; tiny.
- **12b, our own schema code** (§3), optional: after 12a, touching only `src/tools/` and the two disclosure tools;
  parallel with slices 7–9.
