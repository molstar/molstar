# The trace (`src/trace/` of `@molstar/ai`, the `./trace` entry)

A read-only projection of a harness's records into plain data for a trace view (Prism's Inspector: Trace, Tools,
Limits). It shows only what is committed, so it is safe to call while a turn runs; re-read on `client.watch()`, not on
every streamed token. Imports `messages/` and `harness/` only (and may import `session/compose.ts`, which it does not
yet: its markers and label kinds, `LABELS`, are copied).

## Files (`src/trace/`)

| File        | What                                                                                                             |
| ----------- | ---------------------------------------------------------------------------------------------------------------- |
| `read.ts`   | `readTrace`, `TraceProfile`, `createTraceCache`, `defaultRequestOf`, `defaultWindowPart`, `estimatedTokens`      |
| `types.ts`  | the views: `TraceData`, `TurnView`, `StepView`, `WindowView`, `BudgetView`, `TokensView`, tools and limits views |
| `tools.ts`  | `toolsView(input)`: the Tools tab from plain inputs                                                              |
| `format.ts` | `k`, `ms`, `credits`, `cachedShare`, `formatUsage` (`@molstar/ai/trace`)                                         |

## `readTrace(reader: HarnessReader, profile?: TraceProfile): Promise<TraceData>`

Reads `execution()`, `events()` (paged, 500 per page), `turns()`, and the heads, entries and blobs the step events name.
For a model call the request blob is its evidence (`head`, model, `segments`): the head is the window the model saw; the
result blob is its reply. For a tool call the request blob is `{ toolId, arguments }`, the result `{ status, parts }`. A
turn's prompt is the last user message of its begin step's `main` head.

`TraceProfile` (all optional): `requestOf(message)` (prompt and focus; default from the `session:request` label's data,
else the `Request: ` / `Focus: ` markers), `windowPart(message)` (`guidance | tools | state | conversation`; default by
the label for guidance, updates and definitions, else by shape and markers), `toolCall(call)` (title, read, changed:
what the host knows of a call), `credits(model, usage)`, `ceilings` (by budget kind; default the turn's limits),
`toolNames`, `window`, `now`, `cache` (a `TraceCache` kept across re-reads: heads, entries and blobs are immutable, so a
re-read fetches only what is new).

`TraceData`: `harnessId`, `execution` (seq, head, pending step, parked), `budget` (the kind nearest its ceiling in the
latest turn), `preamble` (guidance and tool definitions with token estimates), `turns`, `loaded` and `declared` tools
(from `main`).

- `TurnView`: `id`, `prompt`, `focus?`, `status` (`thinking | working | waiting | done | stopped | failed`), times (ms
  from the first event), `tokens`, `steps`, `milestones`, `failure?`.
- `StepView` is a `ModelStepView` (model, tokens, reply, tool calls, `finish`, `window`, `prefixBreak`) or a
  `ToolStepView` (`toolId`, `name`, `title?`, `read?`, `args`, `status`, `result`, `changed`); each with `attempts` and
  `under?` (the parent call: a nested run, a compaction). Every call of the turn, nested included, in start order.
- `WindowView`: the messages the model saw, each with its part, a token estimate (chars/4) and `fresh` (new since the
  run's previous call).
- `PrefixBreak { reason, kept, of }`: where a request stopped extending the previous one (by provider `segments`): model
  or configuration changed, tools changed, system text changed, or the conversation changed at message N.
- `TokensView`: `input` (uncached + cache read + cache write), `output`, `cached`, `written`, `reasoning`, `estimated`,
  `credits` (display only).

## The other tabs

- `toolsView({ tools: ToolEntry[], skills, loaded, declared, titles })` → `ToolsData`: groups core, each skill (its
  tools, then extras), single; per tool its definition, size in tokens, state (`offered | deferred | loaded`),
  unavailability. The host builds the entries from its toolkit, its placements as `placementOf` gives them (a `group` is
  a skill, an `entry`'s tools show under the entry: pass the entry in `skills` for its title); the trace gives `loaded`
  and `declared`.
- Limits views (`LimitRowView`, `LimitsGroupView`, `LimitsData`) are types only; the host builds them.

## Consumers

`@molstar/prism-assistant`'s Inspector (`src/inspector/`) uses the types and `format`; in the Prism viewer its data
comes from `src/viewer/inspector-sources.ts` (`readTrace` on the session's client, re-read on `watch()`, one read at a
time; `toolsView` from the viewer's toolkit), in the stories from mock data (`src/testing/mock-inspector.ts`).
