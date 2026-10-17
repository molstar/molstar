# The trace (`src/trace/` of `@molstar/ai`, the `./trace` entry)

A read-only projection of a harness's records into plain data for a trace view (Prism's Inspector: Trace, Tools,
Limits). It shows only what is committed, so it is safe to call while a turn runs; re-read on `client.watch()`, not on
every streamed token. Imports `messages/` and `harness/` only: the session's markers and labels (`REQUEST_MARKER`,
`LABEL_REQUEST`, …) are in `messages/markers.ts`, so `@molstar/ai/trace` pulls neither the session nor the tools (the
tools' label kinds in `LABELS` are copied).

## Files (`src/trace/`)

| File        | What                                                                                                             |
| ----------- | ---------------------------------------------------------------------------------------------------------------- |
| `read.ts`   | `readTrace`, `TraceProfile`, `createTraceCache`, `defaultRequestOf`, `defaultWindowPart`, `LABELS`               |
| `types.ts`  | the views: `TraceData`, `TurnView`, `StepView`, `WindowView`, `BudgetView`, `TokensView`, tools and limits views |
| `tools.ts`  | `toolsView(input)`: the Tools tab from plain inputs                                                              |
| `format.ts` | `k`, `ms`, `credits`, `cachedShare`, `formatUsage` (`@molstar/ai/trace`)                                         |

## `readTrace(reader: HarnessReader, profile?: TraceProfile): Promise<TraceData>`

Reads `execution()`, `events()` (paged, 500 per page), `turns()`, and the heads, entries and blobs the step events name.
For a model call the request blob is its evidence (`head`, model, `prefix`): the head is the window the model saw; the
result blob is its reply. For a tool call the request blob is `{ toolId, arguments }`, the result `{ status, parts }`. A
turn's prompt is the last user message of its begin step's `main` head.

`TraceProfile` (all optional): `requestOf(message)` (prompt and focus; default from the `session:request` label's data,
else no request recorded), `windowPart(message)` (`guidance | tools | state | conversation`; default by the label for
guidance, updates, definitions and requests (state when text precedes the request), else by shape), `toolCall(call)`
(title, read, changed: what the host knows of a call), `credits(model, usage)`, `ceilings` (by budget kind; default the
turn's limits), `toolNames`, `window` (the model's context window, for scale), `cache` (a `TraceCache` kept across
re-reads: heads, entries and blobs are immutable, so a re-read fetches only what is new).

`TraceData`: `harnessId`, `execution` (seq, head, pending step, parked), `budget` (the kind nearest its ceiling in the
latest turn), `preamble` (guidance and tool definitions with token estimates), `turns`, `loaded` and `declared` tools
(from `main`).

- `TurnView`: `id`, `prompt`, `focus?`, `status` (`thinking | working | waiting | done | stopped | failed`), times (ms
  from the first event), `tokens`, `steps`, `milestones`, `failure?`.
- `StepView` is a `ModelStepView` (model, tokens, reply, tool calls, `finish`, `window`, `prefixBreak`) or a
  `ToolStepView` (`toolId`, `name`, `title?`, `read?`, `args`, `status`, `result`, `changed`); each with `attempts` and
  `under?` (the parent call: a nested run, a compaction). Every call of the turn, nested included, in start order.
- `WindowView`: the messages the model saw, each with its part, a token estimate (`estimatedTokens` of `messages/`,
  chars/4) and `fresh` (new since the run's previous call).
- `PrefixBreak { reason, kept, of }`: where a request stopped extending the previous one (by the provider's `prefix`):
  model or configuration changed, tools changed, system text changed, or the conversation changed at message N.
- `TokensView`: `input` (uncached + cache read + cache write), `output`, `cached`, `written`, `reasoning`, `estimated`,
  `credits` (display only).

## Incremental reads

A host reads again at every change of the harness while a turn runs, so a read costs what changed (`TraceCache`, kept by
the host across reads of one harness; another harness resets it). Reason: a read that maps every call's window is calls
× messages, and its fresh view objects rebuild every row built on them.

- **Events**: fetched after the last seq held, grouped by turn as they come (one pass over the new ones).
- **A model call** keeps its view, while its `CallRecord` is the same object and it follows the same calls: the run's
  previous call (what is `fresh`), the run's previous request and the latest request of any run (where the provider's
  cache broke). Those are held as the reads of those calls, so a changed call changes every call after it, and nothing
  else. Its window entries are made once per message (`windowEntries`, by entry id). A new call costs its own window.
- **A tool call** keeps its view while its record and the host's `toolCall` notes on it are the same.
- **A turn** keeps its view, and what it carries on to the turns after it (the tool names it learned, the tools first
  loaded in it, the preamble), while its events (count and last), summary and status, the flow into it and the first
  event's time are the same. A turn with a new call is rebuilt from its calls' kept views, and keeps its own view when
  every call is the one it had. So after a turn ends nothing of it is read again, and the running turn costs its calls
  listing plus the new call's window.
- **The definitions in `main`** are read from where the head extended the one read before.
- **Kept only when whole**: a call or turn whose blob or head could not be read is read again, not kept. The host's
  `toolCall` notes of a turn that has not changed are not asked again (a notes change shows with the turn's next event).
- A cold read of N calls is Σ windows: the output itself holds a window per call. The scaling guards
  (`_test/reading.scaling.ts`, `pnpm -C packages/ai test:scaling`) read one more call of N, and N reads with nothing
  new, at n and 4n; a read with no cache is the control.

## The other tabs

- `toolsView({ tools: ToolEntry[], skills, loaded, declared, titles })` → `ToolsData`: groups core, each skill (its
  tools, then extras), single; per tool its definition, size in tokens, state (`offered | deferred | loaded`),
  unavailability. The host builds the entries from its toolkit, its placements as `placementOf` gives them (a `group` is
  a skill, an `entry`'s tools show under the entry: pass the entry in `skills` for its title); the trace gives `loaded`
  and `declared`.
- Limits views (`LimitRowView`, `LimitsGroupView`, `LimitsData`) are types only; the host builds them.

## The Inspector

What the trace is shaped for:

- Three views: **Trace** (an overview strip, the turns as rows, nested runs indented; a model call shows exactly what it
  saw), **Tools** (as the model is offered them: core, skills, single; loaded, unavailable and why; sizes), **Limits**
  (each limit, its value and default; the turn budget against what the latest turn used).
- A page over the whole app, not a panel: it needs the room.
- One usage formatter (`format.ts`) for the Inspector, the Assistant's card and the History.
- Credits are a display unit from the host's prices, never a limit.

`@molstar/prism-assistant`'s Inspector (`src/inspector/`) uses the types and `format`; in the Prism viewer its data
comes from `src/viewer/inspector-sources.ts` (`readTrace` on the session's client, re-read on `watch()`, one read at a
time; `toolsView` from the viewer's toolkit), in the stories from mock data (`src/testing/mock-inspector.ts`).
