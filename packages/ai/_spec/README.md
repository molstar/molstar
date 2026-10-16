# Mol\* AI: the specification

The agent stack under `packages/ai/`: a durable agent loop (the harness), what feeds it (model adapters, tools) and what
a host drives it through (the session). Generic: nothing here imports Mol\*, Prism, React or rxjs. This spec is for
working on or with the code without reading all of it: what each package is, how they bind, the contracts a host relies
on, and where the code is. Why the design is so is in [.v6/ai](../../../.v6/ai/README.md); code comments cite its
sections (`harness.md §6`, `tools-and-session.md §5`), and those numbers are fixed. The code is the truth for how.

| File                         | What it covers                                                                                      |
| ---------------------------- | --------------------------------------------------------------------------------------------------- |
| [harness.md](harness.md)     | the engine: journal, event log, Execution, runner, steps, decisions, the tool handler, storage      |
| [models.md](models.md)       | `Message`, the `ModelAdapter` port, the adapters, the model catalog                                 |
| [tools.md](tools.md)         | `@molstar/ai-tools`: defining tools, field shapes, results, the toolkit (skills), bindings          |
| [session.md](session.md)     | `@molstar/ai-session`: `AgentSession`, the drive's policies, budgets, compaction, `runOnce`, prices |
| [trace.md](trace.md)         | `@molstar/ai-trace`: the trace read from a harness, for the Inspector                               |
| [fibers.md](fibers.md)       | `@molstar/ai-fibers`: durable append-only streams (transcripts), separate from the harness          |
| [narration.md](narration.md) | `@molstar/ai-narration`: text to speech and audio encoding                                          |
| [testing.md](testing.md)     | the conformance suite, the test doubles, how to test a change                                       |

## Packages and layers

```
ai-agent ─┬─ ai-harness ─┬─ ai-components ─┐
          │              ├─ ai-tools ──────┼─ ai-session
          │              ├─ ai-trace       │
          │              └─ ai-conformance │
          └─ ai-model-adapters             (host passes adapters in)
ai-fibers (nothing) ── ai-components (its backends), ai-conformance/fibers
ai-narration (nothing of the workspace)
```

| Package                      | Holds                                                                                                             | Imports of the workspace                                |
| ---------------------------- | ----------------------------------------------------------------------------------------------------------------- | ------------------------------------------------------- |
| `@molstar/ai-agent`          | `Message` and parts, `ToolDefinition`, `ModelConfiguration`, the `ModelAdapter` port, `ModelError`, `hashMessage` | none                                                    |
| `@molstar/ai-harness`        | the contracts (`src/contracts/`), the host protocol, the runner, `defaultDrive`/`drive`, `callTool`               | `ai-agent`                                              |
| `@molstar/ai-components`     | storage (memory; IndexedDB on a durable core); the in-process supervisor; the fibers' IndexedDB rows              | `ai-agent`, `ai-harness`, `ai-fibers`; `node:*` allowed |
| `@molstar/ai-model-adapters` | Anthropic, OpenAI (Responses, chat completions), Google, on the AI SDK; `./model-catalog` (no dependencies)       | `ai-agent`                                              |
| `@molstar/ai-tools`          | `defineTool`, field shapes (`./zod` for zod input), results, bindings, repeated reads, nested tools, the toolkit  | `ai-agent`, `ai-harness`; zod only here                 |
| `@molstar/ai-session`        | `AgentSession`, the session's drive, budgets, composition, `./compaction`, `runOnce`, prices                      | `ai-agent`, `ai-harness`, `ai-components`, `ai-tools`   |
| `@molstar/ai-trace`          | `readTrace` and the views the Inspector draws; `./format`                                                         | `ai-agent`, `ai-harness`                                |
| `@molstar/ai-fibers`         | the fiber kernel and host, its own storage port                                                                   | none                                                    |
| `@molstar/ai-conformance`    | `runHarnessConformance`, `meterStorage`; `./testing` (clock, ids, scripted model and tools); `./fibers`           | `ai-agent`, `ai-harness`, `ai-fibers`                   |
| `@molstar/ai-narration`      | TTS providers (Gemini, OpenAI), WAV encoding; `./mediabunny` for AAC/Opus                                         | none                                                    |

The import rules are Biome's (`noRestrictedImports` per package in the root `biome.json`); they are architecture, so a
failing import means "not here", not "loosen the rule". Each `src/index.ts` is the package's public surface; deep
imports (`./*`) exist but a host should not need them.

## How it binds: one turn, end to end

1. **The host** makes an `AgentSession` (`ai-session`) with a model source (adapters from `ai-model-adapters`, looked up
   by provider), a toolkit (`ai-tools`), its API object, guidance and optionally `HarnessStorage` (`ai-components`;
   memory by default).
2. **`session.turn(request, hooks)`** composes the user message (`[preamble][state] Focus: … notes Request: …`), plus
   guidance and tool definitions on the first turn, and `client.submit`s it as an inbox item whose `options` carry the
   `TurnConfiguration` (model, tool set, budget limits).
3. **The host protocol** (`createHarnessHost(storage)`, `ai-harness`) records `submitted` and, with no runner on the
   harness, returns a **launch key**. The session hands it to the **supervisor** (`createInProcessSupervisor`,
   `ai-components`), which starts a **runner** (`createRunner`) that `claim`s a lease.
4. **The runner** asks the **drive** `next-turn`, `begin`s the turn, then runs **rounds**: a model step (project through
   the adapter, admit against window and budget, stream) and, if the model called tools, a tools step (repair,
   permission, schedule, run handlers, append the assistant message and every result in one move). Each step ends in one
   atomic `commit`. Every choice on the way is a **decision** the drive answers (`createDrive` in the session, over
   `defaultDrive`).
5. **Tool handlers** come from `createToolBindings(toolkit.tools, api, …)`; each tool's `run` gets only the API modules
   it declared. Results are rendered to text for the model; effects (`none`/`applied`/`unknown`) steer retries.
6. **The session** follows `client.watch()` (level-triggered seq) and `runner.observe()` (live deltas) to fire
   `onText`/`onTool`/`onStatus`/`onTurn`, and resolves `TurnResult` when the turn ends or parks (`awaiting-approval`;
   `session.resume({ decisions, hooks })` goes on).
7. **Afterwards** anything can read the harness: `readTrace` (`ai-trace`) for the Inspector, `client.turns()` and
   `events()` for history. A reload with IndexedDB storage and the same `harnessId` finds the conversation; a turn whose
   runner died is taken up after its lease (`session.interrupted()`, `resume()`).

## Where to start, by task

| Task                                       | Read                                                                                                 |
| ------------------------------------------ | ---------------------------------------------------------------------------------------------------- |
| host an agent in an app                    | [session.md](session.md), then [tools.md](tools.md) and [models.md](models.md) § Hosting             |
| add or change a tool                       | [tools.md](tools.md)                                                                                 |
| change how a decision is answered          | [session.md](session.md) § The drive; [harness.md](harness.md) § Decisions                           |
| change the loop, the protocol or a record  | [harness.md](harness.md), then [testing.md](testing.md) (a contract change needs a conformance case) |
| add a provider or a model                  | [models.md](models.md)                                                                               |
| add a storage backend                      | [harness.md](harness.md) § Storage, [testing.md](testing.md)                                         |
| show what happened in a turn               | [trace.md](trace.md)                                                                                 |
| keep a durable, windowed transcript or log | [fibers.md](fibers.md)                                                                               |

## State of the code (2026-10-16)

- Built and tested: everything above. Consumers in the repository: the Prism viewer's agent
  (`@molstar/prism-viewer/agent`: operations as tools, `createViewerAgent` over `AgentSession`, one undo step per turn;
  `.v6/prism-ui/spec/agent.md`) and `@molstar/prism-assistant` (the setup, the turn card, the Inspector over
  `ai-trace`).
- Not built: SQLite backends, more than one active turn, `ai-evals`. What is missing and what was deferred:
  [../TODO.md](../TODO.md); the order: [.v6/ai/plan.md](../../../.v6/ai/plan.md).

## Rules for changes here

- A change to a harness contract comes with a conformance case ([testing.md](testing.md)).
- One home per concern ([.v6/ai/harness.md](../../../.v6/ai/harness.md) §1): if a change moves what a model sees out of
  the journal, evidence out of the event log, or a choice out of the drive, it is the wrong change.
- Tests never call a real model; anything that does costs money, so ask before running it.
- Keep this spec true: a change that moves what a file here says updates it in the same commit, as for `.v6`.
