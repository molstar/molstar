# Mol\* AI

The agent harness and what it needs: durable agent turns, model adapters, storage. Generic and UI-free; Prism hosts it
(`@molstar/prism-assistant` is the Assistant's UI, [../prism-ui/](../prism-ui/)).

## Packages (`packages/ai/`)

| Package                      | Holds                                                                                                                                                     |
| ---------------------------- | --------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `@molstar/ai-agent`          | the model side's shapes: `Message` and its parts, tool definitions, `ModelConfiguration`, the `ModelAdapter` port, `hashMessage`, `digest`                |
| `@molstar/ai-harness`        | the contracts (storage port, supervisor) and the platform-free logic: host protocol, runner, tools, nested runs, drives, observation                      |
| `@molstar/ai-components`     | storage (memory, IndexedDB) on one durable core; the in-process supervisor                                                                                |
| `@molstar/ai-model-adapters` | Anthropic, OpenAI (Responses or chat completions, any compatible endpoint), Google, on the AI SDK; the model catalog (`./model-catalog`, no dependencies) |
| `@molstar/ai-tools`          | tools: input as field shapes, results for the model, bindings, repeated reads, nested tools, schedules, groups (skills) and their disclosure     |
| `@molstar/ai-session`        | an agent session: turns, resume, stop; the default drive; message composition; compaction; one-off runs                                                   |
| `@molstar/ai-conformance`    | the suite every storage backend runs, the storage meter, and test doubles at `./testing` (clock, ids, scripted model and tools)                           |

Opt-in like the Prism modes: nothing in Mol\* or prism-ui imports them, and they import nothing of either.

## Documents

| Document                                           | What it is                                                                                                     |
| -------------------------------------------------- | -------------------------------------------------------------------------------------------------------------- |
| [harness.md](harness.md)                           | the architecture; code comments cite its sections (`harness.md §6`)                                            |
| [packages/ai/\_spec](../../packages/ai/_spec/README.md) | the spec for working with the code: packages, contracts, how they bind, where things are |
| [tools-and-session.md](tools-and-session.md)       | `@molstar/ai-tools` and `@molstar/ai-session`: tools, field shapes, groups, the session, the drive, compaction |
| [plan.md](plan.md)                                 | where it comes from, the move, what is next                                                                    |
| [analysis/viewer-mode.md](analysis/viewer-mode.md) | the agent in the viewer mode: tools, skills, the Inspector, fibers (2026-10-11)                                |

## Decisions

1. **Policy in the drive, mechanism in the harness.** Flow, tool execution, retries and accounting are guaranteed; every
   choice is a decision with a replaceable default, so an app writes only its overrides.
2. **One home per concern**: the journal for what a model sees, the event log for evidence, the drive for choices.
3. **Content-addressed, immutable entries** with no metadata bag: harness bookkeeping never changes a `Message`. What
   an entry _is_ (the request, the guidance, a skill's guidance) is its optional `label`: content, so hashed, never sent
   to a provider, and never read by the harness (2026-10-14).
4. **The runner is trusted.** The host checks leases, compare-and-set, atomicity and replay by key, not the step
   grammar.
5. **Nested runs re-run instead of being durable.** Much smaller than durable child runs; a crashed parent repeats its
   helper.
6. **The watch is level-triggered**, so a dropped connection cannot lose a Stop.
7. **Adapters serve any endpoint of their dialect** under that endpoint's id; replay identity is never a URL.
8. **Platform-free core.** Only `ai-components` touches the platform (IndexedDB); the rest runs in browsers and Node
   alike.
9. **No data migration.** A turn written by another schema version is stopped or reset, not migrated.
