# Mol\* AI

The agent harness and what runs on it: agent turns in memory, model providers, tools, a session. Generic and UI-free;
Prism hosts it (`@molstar/prism-viewer`'s agent, and `@molstar/prism-assistant` for the Assistant's UI;
[../prism-ui/](../prism-ui/)).

## The package (`packages/ai/`)

One package, `@molstar/ai`, in folders under `src/`:

| Folder       | Holds                                                                                                                         |
| ------------ | ----------------------------------------------------------------------------------------------------------------------------- |
| `messages/`  | the shapes both sides share: `Message` and its parts, tool definitions, `ModelConfiguration`, the `Provider` port, `hashMessage`, `digest` |
| `harness/`   | `createHarness`: the in-memory loop (journal, event log, steps, tools, nested runs, drives, observation) and its reader; its cases are `harness/_test/` |
| `providers/` | Anthropic, OpenAI (Responses or chat completions, any compatible endpoint), Google: a thin layer over the AI SDK; the model catalog (`models.ts`) |
| `tools/`     | tools: input as field shapes, results for the model, bindings, repeated reads, nested tools, schedules, the toolkit (skills) and its disclosure |
| `session/`   | `AgentSession`: turns, resume, stop; the session's drive; message composition; compaction (`session/compaction/`)            |
| `trace/`     | the Inspector's trace, a projection of the harness records (with `format.ts`)                                                 |
| `threads/`   | append-only transcripts: the kernel, the host, the row-store port and memory rows; the threads' conformance suite in `_test/` |
| `narration/` | text-to-speech providers behind one port, the audio encoders (WAV; Mediabunny at its own entry)                               |
| `testing/`   | the `createSession` fixture and test doubles (clock, ids, scripted model and tools)                                           |

The entries are the explicit `exports` of `package.json` (no `./*`): `.` (what a host names: messages, the harness's
types, tools, the session with its compaction, threads; nothing that pulls the AI SDK), `./harness`, `./providers`,
`./providers/{anthropic,openai,google}` (one SDK per lazy import), `./models` (the catalog, dependency-free), `./trace`,
`./testing`, `./testing/threads`, `./narration` and `./narration/mediabunny`. The package is `"sideEffects": false`.
The folders keep the layering `messages → harness → tools → session`, with `trace` over the harness and `providers`,
`threads`, `narration` and `testing` beside them; Biome enforces it with per-folder import rules (`biome.json`).

Opt-in: nothing in Mol\* or prism-ui imports it, and it imports nothing of either but `@molstar/core`'s one hash file
(in `messages/`).

## Documents

| Document                                                | What it is                                                                                         |
| ------------------------------------------------------- | -------------------------------------------------------------------------------------------------- |
| [harness.md](harness.md)                                | the architecture; code comments cite its sections (`harness.md §6`)                                |
| [tools-and-session.md](tools-and-session.md)            | `tools/` and `session/`: tools, field shapes, skills, the session, the drive, compaction           |
| [packages/ai/\_spec](../../packages/ai/_spec/README.md) | the spec for working with the code: the folders and entries, how a turn binds them, where things are |
| [plan.md](plan.md)                                      | where it comes from, and what is next                                                              |
| [analysis/viewer-mode.md](analysis/viewer-mode.md)      | the agent in the viewer mode: tools, skills, the Inspector, threads (2026-10-11)                   |
| [simplification/](simplification/README.md)            | the analysis and slices that made the package what it is (2026-10-16)                              |

## Decisions

1. **Policy in the drive, mechanism in the harness.** Flow, tool execution, retries and accounting are guaranteed;
   every choice is a decision with a replaceable default, so an app writes only its overrides.
2. **One home per concern**: the journal for what a model sees, the event log for evidence, the drive for choices.
3. **Content-addressed, immutable entries** with no metadata bag: harness bookkeeping never changes a `Message`. What
   an entry _is_ (the request, the guidance, a skill's guidance) is its optional `label`: content, so hashed, never
   sent to a provider, and never read by the harness. Code reads an entry by its label, never by its text or markers.
4. **One loop, one turn at a time.** The harness is a loop in the page that calls it; each step's writes apply
   together at its end, so a stop, a failure or a reset mid-step leaves the journal as it was. A second prompt waits
   for its own turn: no steering, no interrupts.
5. **Nested runs live in their parent's step**: nothing of them is written before it ends, and a retried parent runs
   them again.
6. **The watch is level-triggered**: a reader re-reads on each seq it hears and misses nothing.
7. **Providers are a thin layer over the AI SDK and serve any endpoint of their dialect** under that endpoint's id;
   replay identity is never a URL. A request is projected once, to `streamText`'s own arguments as JSON (digested and
   sent as it is). The native mid-conversation tool additions stay for the catalog's `appendsTools` models, so the
   prefix stays cacheable as tools are disclosed ([simplification/providers.md](simplification/providers.md)).
8. **Generic and platform-free.** Nothing imports Prism, rxjs or Mol\* (but core's hash file), and nothing touches the
   platform: no storage, no `node:*`. It runs in browsers and Node alike.
9. **In memory only.** No storage port, durability, leases or supervisor: a reload ends the conversation. What an app
   keeps of a conversation is its own (Prism's record is a thread). Persistence was built and taken out: no host
   used it, and with its conformance runs it was a large share of the code and most of the test time
   ([simplification/harness.md](simplification/harness.md)).
10. **Synchronous hashing**: entries, heads and blobs are addressed with `@molstar/core`'s MurmurHash3-128 over UTF-8
    bytes. An id within one page's memory, not a cryptographic hash, so nothing waits on WebCrypto
    ([simplification/harness.md](simplification/harness.md) §7).
11. **Tool inputs are field shapes**, and zod is the package's own: no zod type in the API, no `./zod` entry, and a
    Biome rule restricts zod everywhere outside `packages/ai`
    ([simplification/tool-schemas.md](simplification/tool-schemas.md)).
12. **One package, folder rules.** The folders' layering is architecture, enforced by Biome's `noRestrictedImports` per
    folder: a failing import means "not here", not "loosen the rule".
