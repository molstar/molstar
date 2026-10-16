# Mol\* AI: the specification

The agent stack under `packages/ai/`: an agent loop in memory (the harness), what feeds it (model adapters, tools) and
what a host drives it through (the session). Generic: nothing here imports Prism, React or rxjs, and of Mol\* only
`@molstar/core`'s hash file (in `ai-agent`). This spec is for working on or with the code without reading all of it:
what each package is, how they bind, the contracts a host relies on, and where the code is. Why the design is so is in
[.v6/ai](../../../.v6/ai/README.md); code comments cite its sections (`harness.md §6`, `tools-and-session.md §5`), and
those numbers are fixed. The code is the truth for how.

| File                         | What it covers                                                                             |
| ---------------------------- | ------------------------------------------------------------------------------------------ |
| [harness.md](harness.md)     | the engine: journal, event log, Execution, the loop, steps, decisions, the tool handler    |
| [models.md](models.md)       | `Message`, the `ModelAdapter` port, the adapters, the model catalog                        |
| [tools.md](tools.md)         | `@molstar/ai-tools`: defining tools, field shapes, results, the toolkit (skills), bindings |
| [session.md](session.md)     | `@molstar/ai-session`: `AgentSession`, the drive's policies, budgets, compaction           |
| [trace.md](trace.md)         | `@molstar/ai-trace`: the trace read from a harness, for the Inspector                      |
| [fibers.md](fibers.md)       | `@molstar/ai-fibers`: durable append-only streams (transcripts), separate from the harness |
| [narration.md](narration.md) | `@molstar/ai-narration`: text to speech and audio encoding                                 |
| [testing.md](testing.md)     | the conformance suite, the test doubles, how to test a change                              |

## Packages and layers

```
ai-agent ─┬─ ai-harness ─┬─ ai-tools ──── ai-session
          │              ├─ ai-trace
          │              └─ ai-conformance
          └─ ai-model-adapters          (host passes adapters in)
ai-fibers (nothing) ── ai-conformance/fibers
ai-narration (nothing of the workspace)
```

| Package                      | Holds                                                                                                             | Imports of the workspace                |
| ---------------------------- | ----------------------------------------------------------------------------------------------------------------- | --------------------------------------- |
| `@molstar/ai-agent`          | `Message` and parts, `ToolDefinition`, `ModelConfiguration`, the `ModelAdapter` port, `ModelError`, `hashMessage` | `@molstar/core` (its hash file only)    |
| `@molstar/ai-harness`        | the contracts (`src/contracts/`), `createHarness` (the in-memory loop), `defaultDrive`/`drive`                    | `ai-agent`                              |
| `@molstar/ai-model-adapters` | Anthropic, OpenAI (Responses, chat completions), Google, on the AI SDK; `./model-catalog` (no dependencies)       | `ai-agent`                              |
| `@molstar/ai-tools`          | `defineTool`, field shapes (`./zod` for zod input), results, bindings, repeated reads, nested tools, the toolkit  | `ai-agent`, `ai-harness`; zod only here |
| `@molstar/ai-session`        | `AgentSession`, the session's drive, budgets, composition, `./compaction`                                         | `ai-agent`, `ai-harness`, `ai-tools`    |
| `@molstar/ai-trace`          | `readTrace` and the views the Inspector draws; `./format`                                                         | `ai-agent`, `ai-harness`                |
| `@molstar/ai-fibers`         | the fiber kernel and host, its own storage port                                                                   | none                                    |
| `@molstar/ai-conformance`    | `runHarnessConformance`, `createSession`; `./testing` (clock, ids, scripted model and tools); `./fibers`          | `ai-agent`, `ai-harness`, `ai-fibers`   |
| `@molstar/ai-narration`      | TTS providers (Gemini, OpenAI), WAV encoding; `./mediabunny` for AAC/Opus                                         | none                                    |

The import rules are Biome's (`noRestrictedImports` per package in the root `biome.json`); they are architecture, so a
failing import means "not here", not "loosen the rule". Each `src/index.ts` is the package's public surface; deep
imports (`./*`) exist but a host should not need them.

## How it binds: one turn, end to end

1. **The host** makes an `AgentSession` (`ai-session`) with a model source (adapters from `ai-model-adapters`, looked up
   by provider), a toolkit (`ai-tools`), its API object and guidance.
2. **`session.turn(request, hooks)`** composes the user message (`[preamble][state] Focus: … notes Request: …`), plus
   guidance and tool definitions on the first turn, and begins a turn of its **harness** (`createHarness`, `ai-harness`,
   in memory) with the `TurnConfiguration` (model, tool set, budget limits).
3. **The harness** runs the turn's **begin step** (the drive's `turn-start`: compact first or not; then the messages are
   appended), then **rounds**: a model step (project through the adapter, admit against window and budget, stream) and,
   if the model called tools, a tools step (repair, permission, schedule, run handlers, append the assistant message and
   every result in one move). Each step's writes apply at its end. Every choice on the way is a **decision** the drive
   answers (`createDrive` in the session, over `defaultDrive`).
4. **Tool handlers** come from `createToolBindings(toolkit.tools, api, …)`; each tool's `run` gets only the API modules
   it declared. Results are rendered to text for the model; effects (`none`/`applied`/`unknown`) steer retries.
5. **The session** follows `harness.observe()` (live deltas, each event) to fire `onText`/`onTool`/`onStatus`/`onTurn`,
   and resolves `TurnResult` when the turn's `done` settles: ended, or parked (`awaiting-approval`;
   `session.resume({ decisions, hooks })` goes on).
6. **Afterwards** anything can read the harness through `session.client()` (a `HarnessReader`): `readTrace` (`ai-trace`)
   for the Inspector, `turns()` and `events()` for history. Nothing outlives the page: a reload ends the conversation.

## Where to start, by task

| Task                                       | Read                                                                                                 |
| ------------------------------------------ | ---------------------------------------------------------------------------------------------------- |
| host an agent in an app                    | [session.md](session.md), then [tools.md](tools.md) and [models.md](models.md) § Hosting             |
| add or change a tool                       | [tools.md](tools.md)                                                                                 |
| change how a decision is answered          | [session.md](session.md) § The drive; [harness.md](harness.md) § Decisions                           |
| change the loop or a record                | [harness.md](harness.md), then [testing.md](testing.md) (a contract change needs a conformance case) |
| add a provider or a model                  | [models.md](models.md)                                                                               |
| show what happened in a turn               | [trace.md](trace.md)                                                                                 |
| keep a durable, windowed transcript or log | [fibers.md](fibers.md)                                                                               |

## State of the code (2026-10-16)

- Built and tested: everything above. Consumers in the repository: the Prism viewer's agent
  (`@molstar/prism-viewer/agent`: operations as tools, `createViewerAgent` over `AgentSession`, one undo step per turn;
  `.v6/prism-ui/spec/agent.md`) and `@molstar/prism-assistant` (the setup, the turn card, the Inspector over
  `ai-trace`).
- Not built: persistence of a conversation (removed 2026-10-16 with the storage backends), more than one active turn,
  `ai-evals`. What is missing and what was deferred: [../TODO.md](../TODO.md); the order:
  [.v6/ai/plan.md](../../../.v6/ai/plan.md).

## Rules for changes here

- A change to a harness contract comes with a conformance case ([testing.md](testing.md)).
- One home per concern ([.v6/ai/harness.md](../../../.v6/ai/harness.md) §1): if a change moves what a model sees out of
  the journal, evidence out of the event log, or a choice out of the drive, it is the wrong change.
- Tests never call a real model; anything that does costs money, so ask before running it.
- Keep this spec true: a change that moves what a file here says updates it in the same commit, as for `.v6`.
