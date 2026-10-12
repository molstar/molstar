# Tools and the session

Two generic packages over the harness, ported from Stories Untold's agent (`lib/stories/stories-agent`) with its story
cut out ([analysis/viewer-mode.md](analysis/viewer-mode.md) §4). Code comments cite this document's sections
(`tools-and-session.md §2`).

## 1. `@molstar/ai-tools`: a tool and its result

- **A tool** is
  `defineTool({ id, description, input, modules, optional?, read?, idempotent?, approval?, approvedBy?, run })`. `Api`
  is the host's object of modules (Stories' `StoryToolApi`, Prism's command context): a tool names the modules it uses
  and `run` gets only those, frozen. `run({ args, api, signal, progress, nest?, call })` gets its arguments normalized
  and checked, `call` (`callId`, `key`, `prior`) for idempotency, and returns a `ToolResult`.
- **Input** is an object shape (`ObjectShape`, §2) or, through `@molstar/ai-tools/zod`, a zod object. Either way
  `ai-tools` makes from it, once, the JSON Schema the model sees and the validator; zod is its own business.
- **`ToolResult`** is `{ ok: true, did, changed?, view?, warnings?, body?, next?, contributions?, effect? }` or
  `{ ok: false, code, message, hint?, transient?, effect? }`, rendered for the model always in that order: what it did
  ("…: already so." when nothing had to change), `Changed: …`, `View: …`,
  `Warnings (fix if it affects the request, else tell the user): …`, the body, `Next: …`. Arguments it does not take are
  dropped with a warning naming the ones it takes; a check that fails is one message (the field's path, what it expects,
  what came, its description).
- **Effect** goes on the result frame: `none` (a read, or a failure before any change), `applied`, or `unknown` (a throw
  from a tool that is not a read). A tool may say it; otherwise it is derived as above.
- **After a write**, the bindings' `afterWrite` hook runs before the result is recorded (Stories saves the story there);
  when it throws, the call fails transient with effect `applied`.
- **Repeated reads**: the same read with the same arguments while the host's `version()` is unchanged answers in one
  line ("Unchanged: … act on it"), and six reads in a row get a note to act. Reset per turn.
- **Nested tools**: a tool runs a helper agent under its own call (`nest`), the helper's tools seen by nothing else.
- **Schedule**: `toolSchedule(calls, reads)`: consecutive reads one stage, every other call its own, in the model's
  order.
- **Direct run**: a tool called from code, checked the same way, its result as data.

## 2. Field shapes

```ts
type FieldShape = (
  | { kind: 'string'; enum?: readonly string[]; verbatim?: boolean; format?: string }
  | { kind: 'number' | 'integer'; min?: number; max?: number }
  | { kind: 'boolean' }
  | { kind: 'array'; item: FieldShape; max?: number }
  | ObjectShape
  | { kind: 'variant'; tag: string; variants: Readonly<Record<string, ObjectShape>> }
  | { kind: 'json' }
) & { optional?: boolean; description?: string; default?: Json };
interface ObjectShape { kind: 'object'; fields: Readonly<Record<string, FieldShape>> }
```

Small models send sloppy arguments; before the check, against the shape: `""` and `null` are absent; strings are trimmed
(not `verbatim` ones); `"3"` is a number and `"true"` a boolean where one is asked for, a number a string where a string
is; a list of strings is one comma-joined string where one string is asked for, and `"A, B"` a list where a list of
non-objects is; an enum matches in any case; objects in a list are read by their own fields. A zod object is read as a
field shape for this. Prism maps PD params to field shapes.

## 3. Tool groups (skills)

A **toolkit** holds the registrations, the **core** (loaded from the first turn) and the **groups**: a group (a skill)
is `{ id, title, summary, tags?, guidance, tools, extras?, unavailable?(api) }`; its tools load with its guidance, its
extras by name or when called. Every tool has one place (core, a group's tools or extras, or single), which the toolkit
checks. `disclose_skills` and `disclose_tools` (names configurable) load groups and single tools as contributions: the
definitions as one system message of `tool-definition` parts, the guidance as a user message keyed by the group and a
fingerprint of its text (a changed text is given again). Their descriptions list what there is, generated from the
toolkit. Other disclosable texts (Stories' data sources and servers) are entries `{ id, summary, text, tools? }`.

## 4. `@molstar/ai-session`: a conversation

`AgentSession` is one conversation over one harness: host, client, a runner per launch under an in-process supervisor,
storage (memory by default; IndexedDB to keep it; a host may bring its own `HarnessStorage`).

- `turn(request, hooks)` with `request = { prompt, state?, focus?, notes?, attachments?, model?, preamble? }` →
  `TurnResult` (`completed`, `stopped`, `failed`, `awaiting-approval` with the calls; text, tool events, usage, context
  size, compactions). One turn at a time. `resume(hooks, decisions)`, `stop()`, `interrupted()`, `forget()` (back to an
  empty conversation), `dispose()`.
- **Messages**: the first turn sends the guidance, then the tool definitions (core loaded, the rest declared), then the
  request; later turns send what changed since (new or changed definitions, changed guidance as a user message), then
  the request. The request is `[preamble][state]\n\nFocus: [focus]\n\n[notes]Request: [prompt]`; the markers are
  exported for compaction and the trace.
- **Hooks**: `onText`, `onTool` (started, progress, succeeded, failed, skipped; with `read` and times), `onStatus`
  (`{ kind: model | tools | permission | retrying | widened | budget-extended | compacting, detail? }`), `onTurn`
  (begun, with the turn id, before the first model call; or resumed: where a host takes its checkpoint).
- **Models**: an adapter per provider (a lookup), the default model, windows and output limits the app gives for models
  an adapter does not know, thinking.

## 5. The drive

Defaults over `defaultDrive`, each replaceable (`drive` option):

| Decision     | Policy                                                                                                                                                                                                                                                                                                                                                             |
| ------------ | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ |
| permission   | an answer given on resume, else `approve(call)`, else the tool's `approvedBy`, else park; a denial is told to the model ("The person did not allow this call.")                                                                                                                                                                                                    |
| schedule     | `toolSchedule`                                                                                                                                                                                                                                                                                                                                                     |
| failure      | tool: retry a transient failure (2 times, linear back-off) when its effect is `none` or the tool is idempotent; a failed read is reported; a failed write skips the calls after it in that response. Model: retry a transient failure (3 times) unless the provider asks to wait over 60 s; else fail in plain words. Projection: compact on `capacity`, else fail |
| model-result | cut short at the output limit: once per turn, ask again with twice the room (within the model's maximum)                                                                                                                                                                                                                                                           |
| budget       | extend the short kind by its initial amount toward its ceiling; a nested run fails                                                                                                                                                                                                                                                                                 |
| repair       | not loaded: answer failed, load its definition, skip the rest; unknown: point to the disclosure tools; not enabled, no handler: say so, load nothing                                                                                                                                                                                                               |
| boundary     | compact at the trigger, else proceed (never steer: each prompt is its own turn)                                                                                                                                                                                                                                                                                    |
| next-turn    | begin with the first waiting item, compacting first from the soft threshold                                                                                                                                                                                                                                                                                        |

**Budget**: tokens by kind, the same for every model, scaled with the compaction trigger (80% of the window): tool calls
60 → 200, uncached input max(500k, 4×trigger) → max(2M, 20×), cache reads max(5M, 40×) → max(25M, 200×), output 200k →
1M. **Prices** are the app's (a table per model); credits are a display unit (1 credit = 1,000 uncached input tokens on
a reference model), never a limit.

## 6. Compaction (`@molstar/ai-session/compaction`)

One nested model call over the whole context asks what to carry over; the new head is the guidance, the kept tools, the
carry-over, the memory, the kept groups' guidance, the recent turns as they were and the current turn. It runs at the
trigger mid-turn, from the soft threshold at a turn's start, and forced (a mechanical digest, no model call) when the
context no longer fits; not again until the context has grown by the headroom. A wish carried over must quote the
person's own words. What it keeps is the host's (`CompactionDomain`: what a request is, the core, the groups, the
memory, pins).

## 7. One-off runs

`runOnce({ models, model, guidance, toolkit, api, settlement, prompt, signal, timeoutMs })`: a throw-away session with
terminal tools (`terminalToolFor<Api>()(settlement, { id, input, modules, run })`) whose first call settles the run with
its value and stops the session (the command menu's assistant proposing a filled command). It returns
`{ status: settled | completed | failed | timeout | aborted, value?, text, error? }` and never throws.
