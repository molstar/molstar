# The session (`src/session/` of `@molstar/ai`)

One conversation over one harness, with the policies an app would otherwise write: composition of requests, the drive,
budgets, compaction, approvals, resume. This is what a host uses; it wires everything in
[README.md § How it binds](README.md#how-it-binds-one-turn-end-to-end). Design:
[.v6/ai/tools-and-session.md](../../../.v6/ai/tools-and-session.md) §4–7.

## Files (`src/session/`)

| File          | What                                                                                                        |
| ------------- | ----------------------------------------------------------------------------------------------------------- |
| `session.ts`  | `AgentSession`: the harness (made on first use), `turn`, `resume`, `stop`, hooks from observations, results |
| `types.ts`    | `TurnRequest`, `TurnHooks`, `TurnResult`, `ToolEvent`, `StatusEvent`, `TurnEvent`, `ModelChoice`            |
| `compose.ts`  | `turnInput`, `requestOf` and the markers (`REQUEST_MARKER`, `FOCUS_MARKER`, …)                              |
| `drive.ts`    | `createDrive` (the session's decisions), `turnBudget`, `extendBudget`, `plainFailure`, `DENIED`             |
| `compaction/` | `runCompaction`, `planCompaction`, `compactionFor`, `CompactionDomain`, `toolkitDomain` (`./compaction`)    |

## Making one

```ts
const session = new AgentSession<MyApi>({
  models: { adapterFor: (provider) => adapters[provider], default: { provider: 'anthropic', model: 'claude-sonnet-5-5' } },
  toolkit,                       // createToolkit(...) from tools/
  api,                           // the host's modules; the session adds `disclosure` if absent
  guidance,                      // the system text of the first turn
});
```

Other options (`AgentSessionOptions`, `session.ts`): `startTools` (`'core'` default, or `'all'`), `clock`,
`contextWindow`, `compaction` (overrides or `false`), `compactionDomain`, `maxToolCalls` (200), `thinking` (one
`ReasoningParameters` or per model; pass the catalog's `defaultReasoning`), `maxOutputTokens` (16384), `version()`
(enables repeated-read answers), `afterWrite` (where a host saves after a write), `price`, `drive` (the host's own
`DriveHandlers`, over the session's), and the `DrivePolicy` fields: `approve(call)`, `retries` (3), `retryDelayMs`
(1500), `toolRetries` (2), `toolRetryDelayMs` (1000), `providerNames`. Also `notesHeading` (heads the notes; default
`NOTES_HEADING`).

The harness (made on first use) is `createHarness({ models, tools, drive, clock, price })`: `ModelBindings` over
`adapterFor`, `createToolBindings(toolkit, api, { repeats, nested: NestedTools, afterWrite })` (resolved at call time)
wrapped to observe calls, and `createDrive(...)`; the session subscribes to its observations. It lives in memory, as
long as the session (a reload ends the conversation).

## Turns

- `turn(request, hooks)` → `Promise<TurnResult>`; one at a time (a second throws). `TurnRequest`: `prompt`, `state?`,
  `focus?`, `notes?`, `preamble?`, `attachments?` (`Reference`s with their bytes), `model?` (a `ModelChoice`
  `{ provider, model, thinking?, contextWindow?, maxOutput? }`).
- **Composition** (`compose.ts`): the user message is
  `[preamble]\n\n[state]\n\nFocus: [focus]\n\nNotes since your last turn:\n- [note]…\n\nRequest: [prompt]` (absent parts
  dropped), plus the attachments. The first turn sends before it the guidance (a system message) and the tool
  definitions (core loaded, the rest declared deferred; one system message). Later turns send only what changed since
  (`catchUp`): new or changed definitions as a system message, changed guidance as a user message headed
  `GUIDANCE_UPDATE`.
- **Labels**: every message the session writes is labelled (`LABEL_REQUEST` with `{ prompt, focus? }` as data,
  `LABEL_GUIDANCE`/`LABEL_GUIDANCE_UPDATE` keyed by the text's fingerprint, `guidanceLabel(text)`, and the toolkit's
  `tools:definitions`). Compaction, `catchUp` and the trace read labels; the markers stay in the text for the model and
  are read only for entries written before labels.
- The turn begins with `harness.turn({ messages, configuration })` and the session awaits its `done`. Tool set: the
  toolkit's `toolSet()` as the turn begins (every registration enabled, all but `approval` ones preauthorized, and the
  tools lazy skills will bring; tools.md, "Dynamic").
- **`TurnResult`**: `status` (`completed`, `stopped`, `failed`, `reset`, `awaiting-approval`), `approval?` (pending
  calls), `text` (the last reply), `turnText` (every reply of the turn: one per model call, a retry counted once, a
  widened reply replacing the cut one; for a resumed turn, since the resume), `error?`, `tools: ToolEvent[]`, `usage?`
  (`input`, `output`, `cached`, `written`, `modelCalls`), `context?` (tokens after the turn), `compactions?`.
- **Hooks** (`TurnHooks`): `onText(textSoFar)`, `onTool(ToolEvent)` (top-level calls: `started`, `progress`,
  `succeeded`, `failed`, `skipped`; with `read`, times),
  `onStatus({ kind: model | tools | permission | budget | retrying | widened | budget-extended | compacting, detail? })`,
  `onTurn({ kind: 'begun', turnId, resumed? })` (before the first model call: where a host takes its undo checkpoint).
- **Approval**: a turn parks at `permission`; the result is `awaiting-approval`.
  `resume({ decisions: { [callId]: 'allow' | 'deny' }, hooks? })` continues it. A denial is told to the model as
  `DENIED` ("The person did not allow this call.").
- `stop()` (the running or parked turn ends stopped; changes already applied stay; a stop while the turn is being
  composed keeps it from beginning), `forget()` (reset `main` to the empty head: a fresh conversation, history kept in
  the log; a parked turn ends reset), `dispose()` (a turn still running or parked ends stopped), `client()` (the
  harness's `HarnessReader`, for a trace and the Inspector), `isFresh()`, `contextWindow()`, `budget(choice?)`.
- Gone 2026-10-16: the options `storage`, `harnessId`, `leaseMs`, `interrupted()` and the takeover of a turn after a
  reload (nothing outlives the page).

## The drive (`createDrive`, `drive.ts`)

`drive(hostOverrides, drive(sessionPolicies, defaultDrive()))`: a host handler answering `'pass'` leaves the decision to
the session's.

| Decision       | The session's answer                                                                                                                                                                                                                                                                                                                                                                                                                            |
| -------------- | ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `permission`   | the decision given on resume, else `approve(call)`, else the tool's `approvedBy`, else `park`; `deny` → `{ deny: DENIED }`                                                                                                                                                                                                                                                                                                                      |
| `schedule`     | `toolSchedule`: a run of consecutive reads together, every other call alone, in the model's order                                                                                                                                                                                                                                                                                                                                               |
| `failure`      | tool: retry transient ones (`toolRetries`, linear back-off) when the effect is `none` or the tool is idempotent; a failed read is reported; a failed write also skips the calls after it. Model: retry transient ones (`retries`, not before the provider's `retryAfterMs`) unless it asks to wait over 60 s (its `retryAfterMs`, else its message); else fail with `plainFailure` words. Projection: `capacity` → forced compaction, else fail |
| `model-result` | cut at the output limit: once per turn, `replace` with twice `maxOutputTokens` (within the model's max); nested runs accept                                                                                                                                                                                                                                                                                                                     |
| `budget`       | `extendBudget`: the short kind grows by its initial amount toward its ceiling; nested runs fail                                                                                                                                                                                                                                                                                                                                                 |
| `repair`       | not loaded: fail it, load its definition (a contribution), skip the rest of the response; unknown (or not loaded and not registered): point to the disclosure tools; not enabled or no handler: "cannot run now"; invalid arguments: "takes its arguments as one JSON object"                                                                                                                                                                   |
| `boundary`     | compact at the trigger, else proceed                                                                                                                                                                                                                                                                                                                                                                                                            |
| `turn-start`   | compact first from the soft threshold, else proceed                                                                                                                                                                                                                                                                                                                                                                                             |

**Budget** (`turnBudget(trigger, maxToolCalls)`): tokens by kind, the same for every model, scaled with the compaction
trigger. Initial → ceiling: tool calls 60 → `maxToolCalls` (200), uncached input max(500k, 4×trigger) → max(2M,
20×trigger), cache reads max(5M, 40×trigger) → max(25M, 200×trigger), output 200k → 1M.

## Compaction (`src/session/compaction/`, also from the root entry)

- **When**: `trigger` (mid-turn at `boundary`, the request reaching `triggerInputTokens`), `soft` (a turn's start, from
  `softShare` of the trigger), `force` (a `capacity` failure: a mechanical digest, no model call). Not again until the
  context has grown by the headroom since the last one (`lastCompaction` from the event log).
- **Sizes** (`compactionFor(window)`): trigger 80% of the window, soft 85% of the trigger, target 60% of the soft
  threshold, carry-over ≤ 8k tokens, recent turns kept verbatim up to a tenth of the trigger (≤ 64k).
- **How** (`runCompaction(context, options)`): one nested model call (`context.run`, no tools) over the context asks
  what to carry over (`CARRY_OVER_HEADER`; a carried wish must quote the person's words, `checkWishes`);
  `planCompaction` builds the `HeadEdit`: the latest guidance (a later update included), kept tools, the carry-over, the
  memory, kept groups' guidance, recent turns as they were, the current turn; each new message labelled
  (`session:guidance`, `tools:definitions`, `LABEL_CARRY_OVER`, `LABEL_MEMORY`). If the call fails, a digest stands in.
  Sizes (`sizeOf`, `estimatedTokens`) count what is sent, not labels.
- **What it keeps is the host's**: `CompactionDomain` (`isRequest`, `requestText?`, `coreToolIds`, `groups`, `groupsIn`,
  `groupEntry`, `memory?`, `pins?`, `wishSources?`); the session builds it from the toolkit (`toolkitDomain`) and merges
  `compactionDomain`.

## Gotchas

- A session is one conversation: one harness, in memory.
- Model choices are the host's: the session knows no models; context windows and output limits come from the choice, the
  adapter's capabilities, or `contextWindow`/`maxOutput` given for custom models.
- `onText` streams one model call at a time (it restarts with each top-level call); `TurnResult.text` is the last
  non-empty text, else the last assistant message; `turnText` is the whole turn's.
- Only "a turn is already running" throws from `turn()`; a missing adapter or unknown model is `status: 'failed'`.
  `resume()` returns `undefined` when busy or with nothing to resume; `forget()` throws while a turn runs.
- Approval decisions passed to `resume` are kept for the conversation (a parked turn resumed without one for a call
  parks again). The pending calls are in `TurnResult.approval`.
- Prompts during a running turn are not steered in; the host queues them (Stories kept one follow-up).
- `TurnResult` and hooks cover top-level calls; a helper's calls (nested runs, compaction) are only in the trace.
- The session keeps no person-facing history: a durable transcript is the application's to compose (from the hooks and
  `TurnResult`, into a thread, [threads.md](threads.md)).
