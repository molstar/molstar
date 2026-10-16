# Plan

## Where it comes from

The harness was rewritten in Stories Untold (`lib/anystar`, branch `harness-next`, built 2026-10-05) against a frozen
spec, conformance first, and Stories ran on it from then on. On 2026-10-11 the five packages moved here unchanged in
logic, as `@molstar/ai-*` (one package, `@molstar/ai`, since 2026-10-16: the names below are the history); Stories' own code (its tools, drive and agent) stayed behind. What carried over from that
work:

- **Conformance first.** A contract change comes with a conformance case, never a quiet edit; the suite runs on every
  backend.
- **One home per concern** ([harness.md](harness.md) §1). A change that moves a concern to another home is the thing to
  look for in review.
- **Import rules fail the build** (`biome.json`): nothing imports Mol\*, Prism or rxjs; `messages/` imports nothing of
  the other folders; `harness/` only `messages/`; no `node:*` anywhere (these were package rules before the nine packages became one, 2026-10-16). The original stack had drifted one
  reasonable local change at a time until concerns had no clear home; a rule that fails says "not here".

## 1. The move (2026-10-11)

Done:

- **The workspace**: `packages/ai/*`, the catalog (the AI SDK pinned to its major, zod as its peer), Node 22.13 or later
  everywhere (raised for `node:sqlite`, which is deferred; the minimum stays), the import rules (`biome.json`),
  prism-ui's opt-in boundary.
- **The code**, logic unchanged: relative imports to `.js`, the header, the repository's formatting and lint, `lib`
  ES2024, Vitest to Jest (the suite's `expect(value, message)` became `expectWith`), the old spec's citations pointed at
  [harness.md](harness.md)'s sections.
- **The tests run with the repository's** (`pnpm test:unit`, one Vitest run, ESM throughout). They ran apart while the
  repository was on Jest (CommonJS, the AI SDK is ESM only); that config is gone.
- **The docs**: this folder, short; the spec, its decisions log and the rewrite's plan stay in Stories Untold's archive.

## 2. Next

In this order (reordered 2026-10-12; [analysis/viewer-mode.md](analysis/viewer-mode.md) §7):

1. **The tool layer and the session** (done): `@molstar/ai-tools` and `@molstar/ai-session`, ported from Stories' agent
   with its story cut out ([tools-and-session.md](tools-and-session.md)).
2. **`@molstar/ai-trace` (done) and the Inspector** (Trace, Tools, Limits; a page over the whole app, `ShellPage`), on
   mock data until a session feeds them.
3. **`@molstar/ai-fibers`** (done; now `threads/` in `@molstar/ai`): the kernel, the host, the backends in `ai-components`, the conformance.
4. **PD → field shape** (done: `paramsShape`, `@molstar/prism-assistant/tools/params`), for commands as tools.
5. **A model catalog** (done) in `@molstar/ai-model-adapters/model-catalog`: entries by provider, families matched by
   prefix at a boundary and models as themselves and their dated snapshots, a model's fields over its families' (a new
   version of a family gets the limits, reasoning and quirks; only a listed model has a label and a price). Each holds
   its label and note, limits, reasoning (modes, efforts, the default), dated prices, and quirks (tool additions). The
   adapters read limits and quirks; the session takes `thinking` per model from the host (`defaultReasoning`), so it
   knows no models (`thinkingFor` is gone); a host lists the models and prices calls (`dollarsOf`; credits are its own).
   The host adds its own endpoints and models (a catalog of its own) and what the person types. Gemini 3's rules
   (reasoning, media in tool results) stay a version rule in the Google adapter.
6. **Pause**: one plan with the Prism UI's, and the sequence viewer. Then Prism hosts it (the integration step of
   [../prism-ui/plan.md](../prism-ui/plan.md): commands as tools, `execute` as the handler, the drive, origins, one undo
   step per turn through the mode's checkpointer, [harness.md](harness.md) §18), the skills, `ai-toolscript` and
   `mvs-code`, the MCP server over the same tools.

Also:

- **`@molstar/ai-narration`** (done): narration voiced, ported from Stories' `stories-narration` for the story mode:
  text-to-speech providers behind one port (Gemini, OpenAI), the audio encoded behind another (WAV built in; Mediabunny
  with the browser's own AAC, else Opus, at `./mediabunny`). No MP3 encoding (every MP3 encoder is LAME- or
  shine-derived, so LGPL); a provider's MP3 is kept as it arrives. Loaded lazily in the app (a lazy contribution).

Later:

- **Evals**: `@molstar/ai-evals`, the common part of Stories' evals (`evals/`, `apps/eval`: scenarios run against real
  models, metering, baselines, the trace of each run), with each mode's scenarios in its `agent` folder.
- **Vitest** for the whole repository, ESM throughout (`.v6/plans/checklist.md`); done, `jest:ai` folded in.
- **The suite waits on conditions** (done 2026-10-16, simplification slice 11). The cases that counted event-loop
  rounds of real work (timers, WebCrypto digests, fake IndexedDB) went with the storage-cost and liveness cases, the
  IndexedDB run and the asynchronous hashing; the AI tests left the `serial` Vitest project and run in the parallel
  `unit` project with everything else. Measured: 9 whole-repository runs and 13 `vitest run packages/ai` runs on a
  machine at a load average of 30 to 70 (other jobs, plus a full run alongside), all green; the 53 AI files sum to
  about 4.4 s of 505 s of file time (under 1%).
- **SQLite storage** for the harness and threads (servers, headless), removed 2026-10-12 until a host needs it; the
  IndexedDB and memory backends and the conformance suites are the pattern.
