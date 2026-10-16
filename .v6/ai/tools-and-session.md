# Tools and the session

Two generic folders of `@molstar/ai` over the harness, ported from Stories Untold's agent (`lib/stories/stories-agent`) with its story
cut out ([analysis/viewer-mode.md](analysis/viewer-mode.md) §4). Code comments cite this document's sections
(`tools-and-session.md §2`).

## 1. `tools/`: a tool and its result

- **A tool** is
  `defineTool({ id, description, input, modules, optional?, read?, idempotent?, approval?, approvedBy?, run })`. `Api`
  is the host's object of modules (Stories' `StoryToolApi`, Prism's command context): a tool names the modules it uses
  and `run` gets only those, frozen. `run({ args, api, signal, progress, nest?, call })` gets its arguments normalized
  and checked, `call` (`callId`, `key`, `prior`) for idempotency, and returns a `ToolResult`.
- **Input** is an object shape (`ObjectShape`, §2), a field shape: no zod in the public API (the owner's decision,
  2026-10-16). `tools/` makes from it, once, the JSON Schema the model sees and the validator; zod is internal to
  `packages/ai` (`tools/`, `providers/`), and a Biome rule keeps it out of every file outside that folder.
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
non-objects is; an enum matches in any case; objects in a list are read by their own fields. Prism maps PD params to field shapes.

## 3. Tool groups (skills)

A **toolkit** holds the registrations, the **core** (loaded from the first turn) and the **groups**: a group (a skill)
is `{ id, title, summary, tags?, guidance, tools, extras?, unavailable?(api) }`; its tools load with its guidance, its
extras by name or when called. Every tool has one place (core, a group's tools or extras, or single), which the toolkit
checks. `disclose_skills` and `disclose_tools` (names configurable) load groups and single tools as contributions: the
definitions as one system message of `tool-definition` parts, the guidance as a user message labelled with the group's key, its id and a
fingerprint of its text (`tools:skill`, `skill:<id>@<fingerprint>`; a changed text is given again). Their descriptions list what there is, generated from the
toolkit. Other disclosable texts (Stories' data sources and servers) are entries `{ id, summary, text, tools? }`.

## 4. `session/`: a conversation

`AgentSession` is one conversation over one in-memory harness (`createHarness`; harness.md §4–§6). Until 2026-10-16 it
also set up a host, a runner per launch under a supervisor, and storage; that went with the harness's durability.

- `turn(request, hooks)` with `request = { prompt, state?, focus?, notes?, attachments?, model?, preamble? }` →
  `TurnResult` (`completed`, `stopped`, `failed`, `reset`, `awaiting-approval` with the calls; the last text and the
  turn's, tool events, usage, context size, compactions). One turn at a time. `resume({ decisions?, hooks? })` (a parked turn), `stop()`, `forget()` (back to an
  empty conversation), `dispose()`.
- **Messages**: the first turn sends the guidance, then the tool definitions (core loaded, the rest declared), then the
  request; later turns send what changed since (new or changed definitions, changed guidance as a user message), then
  the request. The request is `[preamble][state]\n\nFocus: [focus]\n\n[notes]Request: [prompt]`; the markers stay in the
  text for the model, and every message the session writes is labelled (`session:request` with the prompt and focus
  as data, `session:guidance`, `session:guidance-update`, `tools:definitions`), so compaction, `catchUp` and the trace
  read labels, not text (the markers only for entries written before labels).
- **Hooks**: `onText`, `onTool` (started, progress, succeeded, failed, skipped; with `read` and times), `onStatus`
  (`{ kind: model | tools | permission | retrying | widened | budget-extended | compacting, detail? }`), `onTurn`
  (begun, with the turn id, before the first model call; or resumed: where a host takes its checkpoint).
- **Models**: a provider per provider id (a lookup, `providerFor`), the default model, windows and output limits the app
  gives for models a provider does not know, thinking (one for every model, or per model: the host passes the model catalog's default).

## 5. The drive

Defaults over `defaultDrive`, each replaceable (`drive` option):

| Decision     | Policy                                                                                                                                                                                                                                                                                                                                                             |
| ------------ | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ |
| permission   | an answer given on resume, else `approve(call)`, else the tool's `approvedBy`, else park; a denial is told to the model ("The person did not allow this call.")                                                                                                                                                                                                    |
| schedule     | `toolSchedule`                                                                                                                                                                                                                                                                                                                                                     |
| failure      | tool: retry a transient failure (2 times, linear back-off) when its effect is `none` or the tool is idempotent; a failed read is reported; a failed write skips the calls after it in that response. Model: retry a transient failure (3 times, not before the provider's `retryAfterMs`) unless it asks to wait over 60 s; else fail in plain words. Projection: compact on `capacity`, else fail |
| model-result | cut short at the output limit: once per turn, ask again with twice the room (within the model's maximum)                                                                                                                                                                                                                                                           |
| budget       | extend the short kind by its initial amount toward its ceiling; a nested run fails                                                                                                                                                                                                                                                                                 |
| repair       | not loaded: answer failed, load its definition, skip the rest; unknown: point to the disclosure tools; not enabled, no handler: say so, load nothing                                                                                                                                                                                                               |
| boundary     | compact at the trigger, else proceed                                                                                                                                                                                                                                                                                    |
| turn-start   | compact first from the soft threshold, else proceed |

**Budget**: tokens by kind, the same for every model, scaled with the compaction trigger (80% of the window): tool calls
60 → 200, uncached input max(500k, 4×trigger) → max(2M, 20×), cache reads max(5M, 40×) → max(25M, 200×), output 200k →
1M.

## 6. Compaction (`session/compaction/`)

One nested model call over the whole context asks what to carry over; the new head is the guidance, the kept tools, the
carry-over, the memory, the kept groups' guidance, the recent turns as they were and the current turn. It runs at the
trigger mid-turn, from the soft threshold at a turn's start, and forced (a mechanical digest, no model call) when the
context no longer fits; not again until the context has grown by the headroom. A wish carried over must quote the
person's own words. What it keeps is the host's (`CompactionDomain`: what a request is, the core, the groups, the
memory, pins). The new head's guidance is the latest one (an update given after a release included); it labels what it
writes (`session:guidance`, `tools:definitions`, `session:carry-over`, `session:memory`).

## 7. One-off runs

One-off runs were removed on 2026-10-16.
