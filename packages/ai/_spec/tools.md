# Tools (`@molstar/ai-tools`)

How a host declares tools, what the model reads back, and how they become the harness's `ToolBindings`. Design:
[.v6/ai/tools-and-session.md](../../../.v6/ai/tools-and-session.md) §1–3. zod is used only inside this package.

## Files (`tools/src/`)

| File          | What                                                                                      |
| ------------- | ----------------------------------------------------------------------------------------- |
| `tool.ts`     | `defineTool`, `toolFor<Api>()`, `directRun`, `definitionOf`                               |
| `shape.ts`    | `FieldShape`, `ObjectShape`, `ArgsOf`, `normalize` (argument repair)                      |
| `input.ts`    | `ToolInput` (shape, JSON Schema, `parse`), `shapeInput`; `zod.ts`: `zodInput`, `verbatim` |
| `result.ts`   | `ToolResult`, `render`, `failed`, `toFrame`, `Effect`                                     |
| `registry.ts` | `ToolRegistration`, `registered` (the handler pump), `createToolBindings`, `AfterWrite`   |
| `toolkit.ts`  | `createToolkit`: core, groups (skills), entries, the disclosure tools                     |
| `repeats.ts`  | `RepeatedReads`                                                                           |
| `nest.ts`     | `NestedTools`, `Nest`, `NestRequest`: helper agents under a call                          |
| `schedule.ts` | `toolSchedule`                                                                            |
| `closest.ts`  | `suggest`, `closestHint` ("did you mean")                                                 |

## Defining a tool

```ts
const tool = toolFor<MyApi>(); // MyApi: the host's object of modules
export const setColor = tool({
  id: 'set_color',
  description: 'Colors a selection. …',
  input: { kind: 'object', fields: { target: { kind: 'string', format: 'A:10-20' }, color: { kind: 'string' } } },
  modules: ['viewer'],          // run gets only these, frozen
  read: false,                  // reads: scheduled together, effect none, repeated-read guard
  idempotent: true,             // safe to retry after an attempt that may have changed something
  approval: false,              // true: each call waits for the person (not preauthorized)
  run: async ({ args, api, signal, progress, nest, call }) => ({ ok: true, did: 'Colored A:10-20 red.', changed: ['chain A'] }),
});
```

- `ToolSpec` also takes `optional` (modules it does without) and `approvedBy(args, api)` (approves without asking).
- `run` gets `args` normalized and checked (typed from the input), `api` (the declared modules), `signal`,
  `progress(text)` (for the person, never the model), `nest` (when the session made one) and `call` (`callId`, `key`,
  `prior`: a non-idempotent tool answers a `prior: 'possible'` re-run from its earlier effect).
- Arguments the tool does not take are dropped with a warning naming what it takes; a failed check is one message
  (`tool.arguments`) the model can act on.
- `directRun(registration)` calls a tool from code: it gives `(args, api, signal, call?) → Promise<ToolResult>` (the
  result as data, no progress, no `nest`), or `undefined` for a tool not made by `defineTool`. Arguments are normalized
  and checked as for a model's call, but unknown ones fail (`tool.arguments`) instead of being dropped.
- `registered(definition, modules, execute, options)` registers a tool with a hand-written definition.

## Input: field shapes

`FieldShape`: `string` (`enum?`, `verbatim?`, `format?`), `number`/`integer` (`min?`, `max?`), `boolean`, `array`
(`item`, `max?`), `object` (`fields`), `variant` (`tag`, `variants`), `json`; each with `optional?`, `description?`,
`default?`. A tool's input is an `ObjectShape` or `zodInput(z.object(…))` from `@molstar/ai-tools/zod`. Either way the
JSON Schema (draft 2020-12) and the validator are made once.

`normalize` repairs what small models send, before the check: `""` and `null` are absent; strings trimmed (not
`verbatim`); `"3"` → number and `"true"` → boolean where asked; a number → string where a string is asked; a list of
strings → one comma-joined string, `"A, B"` → a list; enums match in any case. Prism maps PD params to field shapes
(`@molstar/prism-assistant` `src/tools/params.ts`, `paramsShape`).

## Results

`ToolResult` is `{ ok: true, did, changed?, view?, warnings?, body?, next?, contributions?, effect? }` or
`{ ok: false, code, message, hint?, transient?, effect? }`. `render` writes it for the model, always in this order:
`did` (as the tool wrote it; by convention "…: already so." when nothing had to change), `Changed: a, b.`, `View: …`,
`Warnings (fix if it affects the request, else tell the user): …`, `body`, `Next: …`. `failed(code, message, hint?)`
makes a failure.

**Effect** (on the result frame; drives retries): a tool may set it; otherwise a successful read is `none`, a successful
write `applied`, the tool's own `{ ok: false }` `none`, a throw `none` for a read and `unknown` otherwise. A thrown
error with `transient: true` becomes `tool.transient`. After a write, `afterWrite` runs before the result is recorded;
if it throws, the call fails transient with effect `applied`.

## Bindings

`createToolBindings(registrations, api, { repeats?, nested?, afterWrite? })` → the harness's `ToolBindings`: resolves a
definition by tool id to `registration.handler(api, options)`, else to a nested helper's tool. The handler pumps
progress frames while `run` works, then one result frame (`toFrame`).

- **`RepeatedReads({ version, streak = 6, exempt? })`**: the same read with the same arguments while `version()` is
  unchanged answers in one line ("Unchanged: … act on it"); after `streak` reads in a row each read says to act.
  `reset()` per turn (the session does it). `exempt` (neither reads nor actions: they don't break a streak) defaults to
  `disclose_skills` and `disclose_tools`; the session passes its toolkit's names, and a host building its own should
  pass `exempt: [toolkit.names.skills, toolkit.names.tools]`.
- **`NestedTools(configure)`**: one per session, bound into the bindings. A tool's
  `nest({ messages, tools, api, model, limits })` runs a helper agent under its call (`context.run`, harness nested run)
  whose tools are its own, preauthorized and seen by nothing else; `model` is the host's model choice, turned into a
  configuration by `configure`.
- **`toolSchedule(calls, reads)`**: consecutive reads are one stage, every other call its own, in the model's order (the
  session's `schedule` answer).

## The toolkit (skills and disclosure)

`createToolkit({ tools, core, groups, entries?, names? })` → `Toolkit` (`Api` must extend `DisclosureApi`: the
disclosure tools use its `disclosure` module and every group's `modules`). Every tool has exactly one place (checked; it
throws naming the tool): the **core** (loaded from the first turn), a **group**'s `tools` (loaded with its guidance) or
`extras` (named by its guidance; loaded by name or when called), an **entry**'s tools (never alone), or **single**
(loaded by name).

- `ToolGroup<Api>` (a skill): `{ id, title, summary, tags?, guidance, tools, extras?, modules?, unavailable?(api) }`.
- `DisclosableEntry`: `{ id, summary, text, tools? }` (other texts the model may load: a data source's notes).
- The toolkit makes two **disclosure tools**, `disclose_skills({ skills })` and `disclose_tools({ tools })` (names
  configurable), whose descriptions are the generated catalog. They answer with **contributions**: the definitions as
  one system message of `tool-definition` parts (`definitionEntries`, labelled `LABEL_DEFINITIONS`) and the guidance as
  a user message (`groupEntry`, labelled `LABEL_SKILL` or `LABEL_ENTRY` with the group's key). They read what is loaded
  through the API's `disclosure` module (`DisclosureModule { loaded(), skills() }`), which the session provides.
- `groupKey(x)` is `skill:<id>@<fingerprint>` or `entry:<id>@<fingerprint>`; `groupsIn(messages)` finds the loaded ones
  by their labels' keys, and only keys of the groups' current text, so a changed guidance text is given again (an entry
  written before labels: by its exact text).
- `coreDefinitions()`, `allDefinitions()`, `placementOf(id)`, `groupEntry(x)`, `disclosureTools()` serve the session and
  compaction.

## Gotchas

- `read: true` is a promise: the result depends only on the arguments and the host's state. It allows concurrent
  scheduling and the repeated-read short-cut.
- A tool needing approval must not be preauthorized: the session's tool set leaves `approval` tools out of
  `preauthorized`, which is what makes the harness ask `permission`.
- Don't validate in the handler with your own schema: the input does it, after `normalize`.
- A tool that is not idempotent and writes must use `call.key`/`call.prior` to avoid doing it twice after a crash.
