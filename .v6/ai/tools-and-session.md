# Tools and the session

Two generic folders of `@molstar/ai` over the harness ([harness.md](harness.md)): `tools/` turns a host's functions
into tools a model can call, and `session/` runs a conversation over one harness with the policies an app would
otherwise write. They came from Stories Untold's agent with its story cut out
([analysis/viewer-mode.md](analysis/viewer-mode.md) §4). Code comments cite this document's sections
(`tools-and-session.md §2`); the numbering is fixed. The working spec:
[packages/ai/\_spec/tools.md](../../packages/ai/_spec/tools.md) and
[session.md](../../packages/ai/_spec/session.md).

## 1. `tools/`: a tool and its result

- **A tool** is
  `defineTool({ id, description, input, modules, optional?, read?, idempotent?, approval?, approvedBy?, run })`. `Api`
  is the host's object of modules (Prism's viewer API): a tool names the modules it uses and `run` gets only those,
  frozen. `run({ args, api, signal, progress, nest?, call })` gets its arguments normalized
  and checked, `call` (`callId`, `key`, `prior`) for idempotency, and returns a `ToolResult`.
- **Input** is an object shape (`ObjectShape`, §2), a field shape. `tools/` makes from it, once, the JSON Schema the
  model sees and the validator. zod does that work inside `tools/` and appears in no exported type, so a host never
  writes a schema twice and never depends on a schema library; a Biome rule keeps zod out of every file outside
  `packages/ai` (decision 11, [README](README.md)).
- **`ToolResult`** is `{ ok: true, did, changed?, view?, warnings?, body?, next?, contributions?, effect? }` or
  `{ ok: false, code, message, hint?, transient?, effect? }`, rendered for the model always in that order: what it did
  ("…: already so." when nothing had to change), `Changed: …`, `View: …`,
  `Warnings (fix if it affects the request, else tell the user): …`, the body, `Next: …`. Arguments it does not take are
  dropped with a warning naming the ones it takes; a check that fails is one message (the field's path, what it expects,
  what came, its description).
- **Effect** goes on the result frame: `none` (a read, or a failure before any change), `applied`, or `unknown` (a throw
  from a tool that is not a read). A tool may say it; otherwise it is derived as above.
- **After a write**, the bindings' `afterWrite` hook runs before the result is recorded (where a host saves); when it
  throws, the call fails transient with effect `applied`.
- **Repeated reads**: the same read with the same arguments while the host's `version()` is unchanged answers in one
  line ("Unchanged: … act on it"), and six reads in a row get a note to act. Reset per turn.
- **Nested tools**: a tool runs a helper agent under its own call (`nest`), the helper's tools seen by nothing else.
- **Schedule**: `toolSchedule(calls, reads)`: consecutive reads one stage, every other call its own, in the model's
  order.
- **Direct run**: a tool called from code, its arguments checked by the same shape, its result as data; stricter on
  one point: an argument the tool does not take fails the run, where a model's call has it dropped with a warning.

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
non-objects is; an enum matches in any case; objects in a list are read by their own fields. Prism maps PD params to
field shapes (`paramsShape`).

## 3. Tool groups (skills)

A **toolkit** holds the registrations, the **core** (loaded from the first turn) and the **groups**: a group (a skill)
is `{ id, title, summary, tags?, guidance, tools, extras?, modules?, unavailable?(api), load?(), pending? }`; its tools
load with its guidance, its extras by name or when called; `modules` are what `unavailable` reads; a lazy skill has
`load` (it brings the real group) and `pending` (the tools it will bring, enabled from the turn's start). Every tool has
one place (core, a group's tools or extras, an entry's tools, or single), which the toolkit checks. `disclose_skills`
and `disclose_tools` (names configurable) load groups and single tools as contributions: the definitions as one system
message of `tool-definition` parts, the guidance as a user message labelled `tools:skill` with the group's key
(`skill:<id>@<fingerprint>`, the fingerprint of its text, so a changed text is given again). Their descriptions list
what there is, generated from the toolkit. Other disclosable texts (a data source's notes) are entries `{ id, summary,
text, tools? }`. The toolkit is dynamic: parts are added and removed between turns (`add`), and a lazy skill loads its
tools when disclosed.

## 4. `session/`: a conversation

`AgentSession` is one conversation over one in-memory harness (`createHarness`; harness.md §4–§6), made on first use and
kept as long as the session.

- `turn(request, hooks)` with `request = { prompt, state?, focus?, notes?, attachments?, model?, preamble? }` →
  `TurnResult` (`completed`, `stopped`, `failed`, `reset`, `awaiting-approval` with the calls; the last text and the
  turn's, tool events, usage, context size, compactions). One turn at a time. `resume({ decisions?, hooks? })` (a
  parked turn), `stop()`, `forget()` (back to an empty conversation), `dispose()`.
- **Messages**: the first turn sends the guidance, then the tool definitions (core loaded, the rest declared), then the
  request; later turns send what changed since (new or changed definitions, changed guidance as a user message), then
  the request. The request is `[preamble][state]\n\nFocus: [focus]\n\n[notes]Request: [prompt]`; the markers stay in the
  text for the model, and every message the session writes is labelled (`session:request` with the prompt and focus
  as data, `session:guidance`, `session:guidance-update`, `tools:definitions`), so compaction, `catchUp` and the trace
  find entries by their labels, not their text (a found entry's text is read for display or carry-over: the state
  before the request's marker, an updated guidance's body).
- **Hooks**: `onText`, `onTool` (started, progress, succeeded, failed, skipped; with `read` and times), `onStatus`
  (`{ kind: model | tools | permission | budget | retrying | widened | budget-extended | compacting, detail? }`), `onTurn`
  (begun, with the turn id, before the first model call; or resumed: where a host takes its checkpoint).
- **Models**: a provider per provider id (a lookup, `providerFor`) and the default model. The session knows no models:
  a window and an output limit come from the provider's `capabilities()`, and thinking is the host's (one for every
  model, or per model: the host passes the catalog's `defaultReasoning`).

## 5. The drive

Defaults over `defaultDrive`, each replaceable (`drive` option):

| Decision       | Policy                                                                                             |
| -------------- | -------------------------------------------------------------------------------------------------- |
| `permission`   | an answer given on resume, else `approve(call)`, else the tool's `approvedBy`, else park; a denial is told to the model ("The person did not allow this call.") |
| `schedule`     | `toolSchedule`                                                                                     |
| `failure`      | tool: retry a transient failure (2 times, linear back-off) when its effect is `none` or the tool is idempotent; a failed read is reported; a failed write skips the calls after it in that response. Model: retry a transient failure (3 times, not before the provider's `retryAfterMs`) unless it asks to wait over 60 s; else fail in plain words. Projection: compact on `capacity`, else fail |
| `model-result` | cut short at the output limit: once per turn, ask again with twice the room (within the model's maximum) |
| `budget`       | extend the short kind by its initial amount toward its ceiling; a nested run fails                 |
| `repair`       | not loaded: answer failed, load its definition, skip the rest; unknown: point to the disclosure tools; not enabled, no handler: say so, load nothing |
| `boundary`     | compact at the trigger, else proceed                                                               |
| `turn-start`   | compact first from the soft threshold, else proceed                                                |

**Budget**: tokens by kind, the same for every model, scaled with the compaction trigger (80% of the window): tool calls
60 → 200, uncached input max(500k, 4×trigger) → max(2M, 20×), cache reads max(5M, 40×) → max(25M, 200×), output 200k →
1M.

## 6. Compaction (`session/compaction/`)

One nested model call over the whole context asks what to carry over; the new head is the guidance, the kept tools, the
carry-over, the memory, the kept groups' guidance, the recent turns as they were and the current turn. It runs at the
trigger mid-turn, from the soft threshold at a turn's start, and forced (a mechanical digest, no model call) when the
context no longer fits; not again until the context has grown by the headroom. A wish carried over must quote the
person's own words. What it keeps is the host's (`CompactionDomain`: what a request is, the core, the groups, the
memory, pins). The new head's guidance is the latest one (a later update included); it labels what it
writes (`session:guidance`, `tools:definitions`, `session:carry-over`, `session:memory`).
