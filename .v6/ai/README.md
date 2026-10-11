# Mol\* AI

The agent harness and what it needs: durable agent turns, model adapters, storage. Generic and UI-free; Prism hosts it
(`@molstar/prism-assistant` is the Assistant's UI, [../prism-ui/](../prism-ui/)).

## Packages (`packages/ai/`)

| Package                          | Holds                                                                                                                                      |
| -------------------------------- | ------------------------------------------------------------------------------------------------------------------------------------------ |
| `@molstar/ai-agent`              | the model side's shapes: `Message` and its parts, tool definitions, `ModelConfiguration`, the `ModelAdapter` port, `hashMessage`, `digest` |
| `@molstar/ai-harness`            | the contracts (storage port, supervisor) and the platform-free logic: host protocol, runner, tools, nested runs, drives, observation       |
| `@molstar/ai-harness-components` | storage (memory, IndexedDB, SQLite at `./storage/sqlite`) on one durable core; the in-process supervisor                                   |
| `@molstar/ai-model-adapters`     | Anthropic, OpenAI (Responses or chat completions, any compatible endpoint), Google, on the AI SDK                                          |
| `@molstar/ai-conformance`        | the suite every storage backend runs, the storage meter, and test doubles at `./testing` (clock, ids, scripted model and tools)            |

Opt-in like the Prism modes: nothing in Mol\* or prism-ui imports them, and they import nothing of either.

## Documents

| Document                 | What it is                                                          |
| ------------------------ | ------------------------------------------------------------------- |
| [harness.md](harness.md) | the architecture; code comments cite its sections (`harness.md §6`) |
| [plan.md](plan.md)       | where it comes from, the move, what is next                         |

## Decisions

1. **Policy in the drive, mechanism in the harness.** Flow, tool execution, retries and accounting are guaranteed; every
   choice is a decision with a replaceable default, so an app writes only its overrides.
2. **One home per concern**: the journal for what a model sees, the event log for evidence, the drive for choices.
3. **Content-addressed, immutable entries** with no metadata bag: harness bookkeeping never changes a `Message`.
4. **The runner is trusted.** The host checks leases, compare-and-set, atomicity and replay by key, not the step
   grammar.
5. **Nested runs re-run instead of being durable.** Much smaller than durable child runs; a crashed parent repeats its
   helper.
6. **The watch is level-triggered**, so a dropped connection cannot lose a Stop.
7. **Adapters serve any endpoint of their dialect** under that endpoint's id; replay identity is never a URL.
8. **Platform-free core.** Only `ai-harness-components` touches the platform (IndexedDB, `node:sqlite`); the rest runs
   in browsers and Node alike.
9. **No data migration.** A turn written by another schema version is stopped or reset, not migrated.
