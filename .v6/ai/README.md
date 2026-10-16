# Mol\* AI

The agent harness and what it needs: agent turns in memory, model adapters, tools, a session. Generic and UI-free; Prism hosts it
(`@molstar/prism-assistant` is the Assistant's UI, [../prism-ui/](../prism-ui/)).

## Packages (`packages/ai/`)

| Package                      | Holds                                                                                                                                                     |
| ---------------------------- | --------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `@molstar/ai-agent`          | the model side's shapes: `Message` and its parts, tool definitions, `ModelConfiguration`, the `ModelAdapter` port, `hashMessage`, `digest` (core's MurmurHash3) |
| `@molstar/ai-harness`        | `createHarness`: an in-memory loop (journal, event log, steps, tools, nested runs, drives, observation) and its reader                                    |
| `@molstar/ai-model-adapters` | Anthropic, OpenAI (Responses or chat completions, any compatible endpoint), Google, on the AI SDK; the model catalog (`./model-catalog`, no dependencies) |
| `@molstar/ai-tools`          | tools: input as field shapes, results for the model, bindings, repeated reads, nested tools, schedules, groups (skills) and their disclosure |
| `@molstar/ai-session`        | an agent session: turns, resume, stop; the default drive; message composition; compaction |
| `@molstar/ai-conformance`    | the harness's conformance cases, the `createSession` fixture, and test doubles at `./testing` (clock, ids, scripted model and tools) |

Also `@molstar/ai-trace` (the Inspector's trace), `@molstar/ai-fibers` (transcripts) and `@molstar/ai-narration`.
`@molstar/ai-components` (storage backends, the supervisor) is gone (2026-10-16): the harness keeps nothing outside the
page's memory ([simplification/harness.md](simplification/harness.md)).

Opt-in like the Prism modes: nothing in Mol\* or prism-ui imports them, and they import nothing of either (but
`@molstar/core`'s one hash file, in `ai-agent`; simplification/harness.md §7).

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
4. ~~The runner is trusted~~: there is no host and no runner apart from the loop (2026-10-16): the harness is one
   in-memory loop that applies each step's writes at its end.
5. **Nested runs are not durable**: they live in their parent's step, and nothing of them is written before it ends.
6. **The watch is level-triggered**: a reader re-reads on each seq it hears and misses nothing.
7. **Adapters serve any endpoint of their dialect** under that endpoint's id; replay identity is never a URL.
8. **Platform-free.** Nothing touches the platform: no storage, no `node:*`; it runs in browsers and Node alike.
9. ~~No data migration~~: nothing is stored (2026-10-16), so there is nothing to migrate.
10. **In memory only** (2026-10-16): no durability, leases, supervisor or storage port; a reload ends the conversation.
    Hashing is synchronous, `@molstar/core`'s MurmurHash3-128 (simplification/harness.md §3, §7).
