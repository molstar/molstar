# Mol\* AI

The agent harness and what it needs: agent turns in memory, model providers, tools, a session. Generic and UI-free; Prism hosts it
(`@molstar/prism-assistant` is the Assistant's UI, [../prism-ui/](../prism-ui/)).

## The package (`packages/ai/`)

One package, `@molstar/ai`, in folders under `src/` ([simplification/packaging.md](simplification/packaging.md)):

| Folder       | Holds                                                                                                                                                      |
| ------------ | ---------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `messages/`  | the model side's shapes: `Message` and its parts, tool definitions, `ModelConfiguration`, the `Provider` port, `hashMessage`, `digest` (core's MurmurHash3) |
| `harness/`   | `createHarness`: an in-memory loop (journal, event log, steps, tools, nested runs, drives, observation) and its reader; its cases are `harness/_test/`      |
| `providers/` | Anthropic, OpenAI (Responses or chat completions, any compatible endpoint), Google: a thin layer over the AI SDK (a request is `streamText`'s own arguments); the model catalog (`models.ts`, no dependencies) |
| `tools/`     | tools: input as field shapes, results for the model, bindings, repeated reads, nested tools, schedules, groups (skills) and their disclosure               |
| `session/`   | an agent session: turns, resume, stop; the default drive; message composition; compaction (`session/compaction/`)                                          |
| `trace/`     | the Inspector's trace, a projection of the harness records (with `format.ts`)                                                                              |
| `threads/`   | transcripts: the kernel, the host, the memory rows; the threads' conformance suite is `threads/_test/`                                                     |
| `narration/` | text-to-speech providers behind one port, the audio encoders (WAV; Mediabunny at its own entry)                                                            |
| `testing/`   | the `createSession` fixture and test doubles (clock, ids, scripted model and tools)                                                                        |

The entries (`package.json` `exports`; explicit, no `./*`): `.` (messages, the harness's types, tools, the session with
its compaction, threads: what a host names, and nothing that pulls the AI SDK), `./harness`, `./providers`,
`./providers/{anthropic,openai,google}` (one SDK per lazy import), `./models` (the catalog, dependency-free), `./trace`,
`./testing`, `./testing/threads` (the threads' conformance suite), `./narration` and `./narration/mediabunny`. The package
is `"sideEffects": false`. The folders keep the layering `messages → harness → {tools, trace} → session`, with
`providers`, `threads`, `narration` and `testing` beside it; Biome enforces it as per-folder import rules (`biome.json`).

`@molstar/ai-components` (storage backends, the supervisor) is gone (2026-10-16): the harness keeps nothing outside the
page's memory ([simplification/harness.md](simplification/harness.md)).

Opt-in like the Prism modes: nothing in Mol\* or prism-ui imports it, and it imports nothing of either (but
`@molstar/core`'s one hash file, in `messages/`; simplification/harness.md §7).

## Documents

| Document                                           | What it is                                                                                                     |
| -------------------------------------------------- | -------------------------------------------------------------------------------------------------------------- |
| [harness.md](harness.md)                           | the architecture; code comments cite its sections (`harness.md §6`)                                            |
| [packages/ai/\_spec](../../packages/ai/_spec/README.md) | the spec for working with the code: packages, contracts, how they bind, where things are |
| [tools-and-session.md](tools-and-session.md)       | `tools/` and `session/` of `@molstar/ai`: tools, field shapes, groups, the session, the drive, compaction |
| [plan.md](plan.md)                                 | where it comes from, the move, what is next                                                                    |
| [analysis/viewer-mode.md](analysis/viewer-mode.md) | the agent in the viewer mode: tools, skills, the Inspector, threads (2026-10-11)                                |

## Decisions

1. **Policy in the drive, mechanism in the harness.** Flow, tool execution, retries and accounting are guaranteed; every
   choice is a decision with a replaceable default, so an app writes only its overrides.
2. **One home per concern**: the journal for what a model sees, the event log for evidence, the drive for choices.
3. **Content-addressed, immutable entries** with no metadata bag: harness bookkeeping never changes a `Message`. What
   an entry _is_ (the request, the guidance, a skill's guidance) is its optional `label`: content, so hashed, never sent
   to a provider, and never read by the harness (2026-10-14). Labels are the only way an entry is read: with nothing persisted, no entry
   predates them, and nothing reads one by its text or markers (2026-10-16).
4. ~~The runner is trusted~~: there is no host and no runner apart from the loop (2026-10-16): the harness is one
   in-memory loop that applies each step's writes at its end.
5. **Nested runs are not durable**: they live in their parent's step, and nothing of them is written before it ends.
6. **The watch is level-triggered**: a reader re-reads on each seq it hears and misses nothing.
7. **Providers are a thin layer over the AI SDK and serve any endpoint of their dialect** under that endpoint's id;
   replay identity is never a URL. A request is projected once, to `streamText`'s own arguments as JSON (retained,
   digested, sent as it is); the native mid-conversation tool additions stay for the catalog's `appendsTools` models
   (2026-10-16, [simplification/providers.md](simplification/providers.md)).
8. **Platform-free.** Nothing touches the platform: no storage, no `node:*`; it runs in browsers and Node alike.
9. ~~No data migration~~: nothing is stored (2026-10-16), so there is nothing to migrate.
10. **In memory only** (2026-10-16): no durability, leases, supervisor or storage port; a reload ends the conversation.
    Hashing is synchronous, `@molstar/core`'s MurmurHash3-128 (simplification/harness.md §3, §7).
