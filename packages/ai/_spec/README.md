# Mol\* AI: the specification

The agent stack, the one package `@molstar/ai` (`packages/ai/`): an agent loop in memory (the harness), what feeds it
(providers, tools) and what a host drives it through (the session). Generic: nothing here imports Prism, React or rxjs,
and of Mol\* only `@molstar/core`'s hash file (in `messages/`). This spec is for working on or with the code without
reading all of it: what each folder is, how they bind, the contracts a host relies on, where the code is, and the
decisions (each file's "Decisions": the rule, then why). Code comments cite numbered sections (`harness.md §6`,
`tools.md §2`, `session.md §3`). The code is the truth for how.

| File                         | What it covers                                                                                 |
| ---------------------------- | ---------------------------------------------------------------------------------------------- |
| [harness.md](harness.md)     | the loop: journal, event log, Execution, steps, decisions, the tool handler, nested runs       |
| [models.md](models.md)       | `messages/` and `providers/`: `Message`, the `Provider` port, the providers, the model catalog |
| [tools.md](tools.md)         | `tools/`: defining tools, field shapes, results, bindings, the toolkit (skills)                |
| [session.md](session.md)     | `session/`: `AgentSession`, composition, the drive's policies, budgets, compaction             |
| [trace.md](trace.md)         | `trace/`: the trace read from a harness, for the Inspector                                     |
| [threads.md](threads.md)     | `threads/`: append-only transcripts over a row store, apart from the harness                   |
| [narration.md](narration.md) | `narration/`: text to speech and audio encoding                                                |
| [testing.md](testing.md)     | the harness cases, the test doubles (`testing/`), how to test a change                         |

## Folders and layers

```
messages ─┬─ harness ─┬─ tools ──── session
          │           ├─ trace
          │           └─ testing
          ├─ providers          (the host passes providers in)
          └─ threads ── testing/threads (the suite)
narration (nothing of the other folders)
```

| Folder (`src/`) | Holds                                                                                                                                                          | Imports of the other folders                 |
| --------------- | -------------------------------------------------------------------------------------------------------------------------------------------------------------- | -------------------------------------------- |
| `messages/`     | `Message` and parts, `ToolDefinition`, `ModelConfiguration`, the `Provider` port, `ModelError`, `hashMessage`; the session's markers and labels (`markers.ts`) | none (`@molstar/core`'s hash file only)      |
| `harness/`      | the contracts (`contracts/`), `createHarness` (the in-memory loop), `defaultDrive`/`drive`; its cases are `_test/`                                             | `messages`                                   |
| `providers/`    | Anthropic, OpenAI (Responses, chat completions), Google, on the AI SDK; `models.ts` (the catalog, no dependencies)                                             | `messages`; the AI SDK only here             |
| `tools/`        | `defineTool`, field shapes, results, bindings, repeated reads, nested tools, the toolkit                                                                       | `messages`, `harness`; zod only here         |
| `session/`      | `AgentSession`, the session's drive, budgets, composition, `compaction/`                                                                                       | `messages`, `harness`, `tools`               |
| `trace/`        | `readTrace` and the views the Inspector draws; `format.ts`                                                                                                     | `messages`, `harness`                        |
| `threads/`      | the thread kernel and host, its row-store port, the memory rows; the suite is `_test/`                                                                         | `messages` (`canonicalJson`, `Clock`)        |
| `testing/`      | `createSession`; the clock, ids, scripted model and tools; the threads' suite (`testing/threads/`)                                                             | `messages`, `harness` (the suite: `threads`) |
| `narration/`    | TTS providers (Gemini, OpenAI), WAV encoding; `mediabunny.ts` for AAC/Opus                                                                                     | none                                         |

The import rules are Biome's (`noRestrictedImports` per folder in the root `biome.json`, over relative imports): zod is
allowed only in `tools/` and `providers/` (and restricted everywhere outside `packages/ai`), the AI SDK only in
`providers/`. They are architecture, so a failing import means "not here", not "loosen the rule". Each folder's
`index.ts` is its public surface.

## The entries

`package.json` `exports`, explicit (no `./*`, and `./_test/*` closed); the package is `"sideEffects": false`.

| Entry                                                       | Source                                    | Holds                                                                                                                            |
| ----------------------------------------------------------- | ----------------------------------------- | -------------------------------------------------------------------------------------------------------------------------------- |
| `@molstar/ai`                                               | `src/index.ts`                            | what a host names: messages, the harness's types, tools, the session with its compaction, threads; nothing that pulls the AI SDK |
| `@molstar/ai/harness`                                       | `src/harness/index.ts`                    | `createHarness`, `emptyHead`, `defaultDrive`, `drive`, the contracts                                                             |
| `@molstar/ai/providers`                                     | `src/providers/index.ts`                  | the three provider factories                                                                                                     |
| `@molstar/ai/providers/anthropic`, `…/openai`, `…/google`   | `src/providers/<name>/provider.ts`        | one factory each, so a lazy import loads one SDK                                                                                 |
| `@molstar/ai/models`                                        | `src/providers/models.ts`                 | the model catalog (data, no dependencies)                                                                                        |
| `@molstar/ai/trace`                                         | `src/trace/index.ts`                      | `readTrace`, `toolsView`, the views, the formatters                                                                              |
| `@molstar/ai/testing`                                       | `src/testing/index.ts`                    | `createSession`, `scriptedModel`, `scriptedTools`, the simulated clock and ids, message builders                                 |
| `@molstar/ai/testing/threads`                               | `src/testing/threads/index.ts`            | `runThreadConformance`, the threads' suite (it imports Vitest; `./testing` does not)                                             |
| `@molstar/ai/narration`, `@molstar/ai/narration/mediabunny` | `src/narration/index.ts`, `mediabunny.ts` | narration; the Mediabunny encoder apart, so Mediabunny stays out of a bundle that does not ask for it                            |

## How it binds: one turn, end to end

1. **The host** makes an `AgentSession` (`session/`) with a model source (providers from `providers/`, looked up by
   provider id), a toolkit (`tools/`), its API object and guidance.
2. **`session.turn(request, hooks)`** composes the user message (`[preamble][state] Focus: … notes Request: …`), plus
   guidance and tool definitions on the first turn, and begins a turn of its **harness** (`createHarness`, `harness/`,
   in memory) with the `TurnConfiguration` (model, tool set, budget limits).
3. **The harness** runs the turn's **begin step** (the drive's `turn-start`: compact first or not; then the messages are
   appended), then **rounds**: a model step (project through the provider, admit against window and budget, stream) and,
   if the model called tools, a tools step (repair, permission, schedule, run handlers, append the assistant message and
   every result in one move). Each step's writes apply at its end. Every choice on the way is a **decision** the drive
   answers (`createDrive` in the session, over `defaultDrive`).
4. **Tool handlers** come from `createToolBindings(toolkit, api, …)`; each tool's `run` gets only the API modules it
   declared. Results are rendered to text for the model; effects (`none`/`applied`/`unknown`) steer retries.
5. **The session** follows `harness.observe()` (live deltas, each event) to fire `onText`/`onTool`/`onStatus`/`onTurn`,
   and resolves `TurnResult` when the turn's `done` settles: ended, or parked (`awaiting-approval`;
   `session.resume({ decisions, hooks })` goes on).
6. **Afterwards** anything can read the harness through `session.client()` (a `HarnessReader`): `readTrace` (`trace/`)
   for the Inspector, `turns()` and `events()` for history.

## What a host relies on

- **One turn at a time.** `turn()` throws only while a turn runs (or is being composed); everything else (a parked turn
  not yet resumed, a missing provider, an unknown model, a failed call) is a `TurnResult` with `status: 'failed'`.
- **A turn always settles**: `completed`, `stopped`, `failed`, `reset`, or `awaiting-approval` until `resume`, `stop` or
  `forget`. `stop()` keeps what tools already applied.
- **The host owns its models, keys and prices**: the session knows no models; a window and an output limit come from the
  provider's `capabilities()` (the catalog, or `modelLimits` for a custom endpoint).
- **Entries are read by their labels** (`session:request`, `tools:skill`, …), never by text: models.md, § Messages.
- **Nothing outlives the page**: the harness is in memory, and a reload ends the conversation. A host that keeps a
  history composes it from the hooks and `TurnResult` (Prism's record is a thread, [threads.md](threads.md)).
- **Reads are async and level-triggered**: re-read on `client().watch()`; a missed emission loses nothing.

## Decisions

- **Generic and platform-free**: nothing imports Prism, React or rxjs, nor Mol\* but `@molstar/core`'s one hash file
  (`messages/hash.ts`); no storage, no `node:*`. It runs in browsers and Node alike, and Prism hosts it, not the other
  way round.
- **One package, folder rules**: one `@molstar/ai` with explicit entries; the layering is Biome's per folder, so a
  failing import means "not here", not "loosen the rule".
- **Opt-in**: nothing in Mol\* or prism-ui imports it; only an app with the assistant pays for it.
- **Lazy where heavy**: each provider (an AI SDK, zod) and Mediabunny is its own entry, so a startup bundle pulls none.

## Where to start, by task

| Task                              | Read                                                                                     |
| --------------------------------- | ---------------------------------------------------------------------------------------- |
| host an agent in an app           | [session.md](session.md), then [tools.md](tools.md) and [models.md](models.md) § Hosting |
| add or change a tool              | [tools.md](tools.md)                                                                     |
| change how a decision is answered | [session.md](session.md) § The drive; [harness.md](harness.md) § Decisions               |
| change the loop or a record       | [harness.md](harness.md), then [testing.md](testing.md) (a contract change needs a case) |
| add a provider or a model         | [models.md](models.md)                                                                   |
| show what happened in a turn      | [trace.md](trace.md)                                                                     |
| keep a windowed transcript or log | [threads.md](threads.md)                                                                 |

## State of the code

- Built and tested: everything above. Consumers in the repository: the Prism viewer's agent (`@molstar/prism-viewer`'s
  `agent/`: operations as tools, `createViewerAgent` over `AgentSession`, one undo step per turn; its record, a thread,
  in `journal/`; `.v6/prism-ui/spec/agent.md`) and `@molstar/prism-assistant` (the setup, the turn card, the Inspector
  over `@molstar/ai/trace`). Narration has no consumer yet.
- Not built: persistence of a conversation, more than one active turn, evals. What is missing and what was deferred:
  [../TODO.md](../TODO.md).

## Rules for changes here

- A change to a harness contract comes with a case in `src/harness/_test/` ([testing.md](testing.md)).
- One home per concern ([harness.md](harness.md) §1): if a change moves what a model sees out of the journal, evidence
  out of the event log, or a choice out of the drive, it is the wrong change.
- Tests never call a real model; anything that does costs money, so ask before running it.
- Keep this spec true: a change that moves what a file here says updates it in the same commit; a citation names a
  section that exists.
