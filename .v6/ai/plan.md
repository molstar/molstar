# Plan

## Where it comes from

The harness was rewritten in Stories Untold (its `lib/anystar`, branch `harness-next`, built 2026-10-05) against a
frozen spec, conformance first, and Stories ran on it from then on. On 2026-10-11 its packages moved here unchanged in
logic, as nine `@molstar/ai-*` packages; Stories' own code (its tools, drive and agent) stayed behind.

On 2026-10-16 they were simplified into what this folder describes ([simplification/](simplification/README.md)): one
package, `@molstar/ai`; a basic in-memory harness in place of the durable one (no storage, leases or supervisor);
providers as a thin layer over the AI SDK; fibers renamed threads; the conformance suite cut to its 59 essential cases;
about 40% of the lines gone.

What carried over from the rewrite:

- **Conformance first.** A contract change comes with a case in `harness/_test/`, never a quiet edit.
- **One home per concern** ([harness.md](harness.md) §1). A change that moves a concern to another home is the thing to
  look for in review.
- **Import rules fail the build** (`biome.json`, per folder): nothing imports Mol\*, Prism or rxjs; `messages/` imports
  nothing of the other folders; `harness/` only `messages/`; no `node:*` anywhere. The original stack had drifted one
  reasonable local change at a time until concerns had no clear home; a rule that fails says "not here".

## 1. The move (2026-10-11)

Done:

- **The workspace**: the packages, the catalog (the AI SDK pinned to its major), Node 22.13 or later everywhere, the
  import rules (`biome.json`), prism-ui's opt-in boundary.
- **The code**, logic unchanged: relative imports to `.js`, the header, the repository's formatting and lint, `lib`
  ES2024, the old spec's citations pointed at [harness.md](harness.md)'s sections.
- **The tests run with the repository's** (`pnpm test:unit`, one Vitest run, ESM throughout).
- **The docs**: this folder, short; the old spec, its decisions log and the rewrite's plan stay in Stories Untold's
  archive.

## 2. Next

In this order (reordered 2026-10-12; [analysis/viewer-mode.md](analysis/viewer-mode.md) §7):

1. **The tool layer and the session** (done): `tools/` and `session/`, ported from Stories' agent with its story cut out
   ([tools-and-session.md](tools-and-session.md)).
2. **The trace and the Inspector** (done): `trace/`, and the Inspector (Trace, Tools, Limits; a page over the whole
   app) in `@molstar/prism-assistant`, on the session's harness.
3. **Threads** (done): `threads/`, append-only transcripts with their own row-store port and conformance suite; the
   viewer's record is one, its IndexedDB rows in `@molstar/prism-viewer`.
4. **PD → field shape** (done: `paramsShape`, now in `@molstar/prism-viewer`'s `agent/params.ts`), for commands as
   tools.
5. **A model catalog** (done) at `@molstar/ai/models`: entries by provider, families matched by prefix at a boundary and
   models as themselves and their dated snapshots, a model's fields over its families' (a new version of a family gets
   the limits, reasoning and quirks; only a listed model has a label and a price). Each holds its label and note,
   limits, reasoning (modes, efforts, the default), dated prices, and quirks (tool additions). The providers read
   limits and quirks; the session takes `thinking` per model from the host (`defaultReasoning`), so it knows no models;
   a host lists the models and prices calls (`dollarsOf`; credits are its own). The host adds its own endpoints and
   models and what the person types. Gemini 3's rules (reasoning, media in tool results) stay a version rule in the
   Google provider.
6. **Prism hosts it** (built: [../prism-ui/spec/agent.md](../prism-ui/spec/agent.md), [harness.md](harness.md) §18):
   the viewer's operations as tools, the drive, origins, one undo step per turn, the skills. Still to come:
   `toolscript` and `mvs-code` (code against the viewer, MVS scripting) and an MCP server over the same tools.

Also:

- **Narration** (done): `narration/`, ported from Stories' `stories-narration` for the story mode: text-to-speech
  providers behind one port (Gemini, OpenAI), the audio encoded behind another (WAV built in; Mediabunny with the
  browser's own AAC, else Opus, at `./narration/mediabunny`). No MP3 encoding (every MP3 encoder is LAME- or
  shine-derived, so LGPL); a provider's MP3 is kept as it arrives. Loaded lazily in the app; no app uses it yet.
- **The review gate** of the simplification ([simplification/plan.md](simplification/plan.md), slice 13): the branch
  reviewed and measured against the analysis before it lands.

Later:

- **Evals**: the common part of Stories' evals (scenarios run against real models, metering, baselines, the trace of
  each run) as a folder or package here, with each mode's scenarios in its `agent` folder (the viewer's first ones are
  its own, `agent/eval.ts`).
- **Vitest** for the whole repository, ESM throughout (`.v6/plans/checklist.md`); done, `jest:ai` folded in.
- **The suite waits on conditions** (done 2026-10-16, simplification slice 11). The cases that counted event-loop
  rounds of real work (timers, WebCrypto digests, fake IndexedDB) went with the storage-cost and liveness cases, the
  IndexedDB run and the asynchronous hashing; the AI tests left the `serial` Vitest project and run in the parallel
  `unit` project with everything else. Measured: 9 whole-repository runs and 13 `vitest run packages/ai` runs on a
  machine at a load average of 30 to 70 (other jobs, plus a full run alongside), all green; the 53 AI files sum to
  about 4.4 s of 505 s of file time (under 1%).
- **Persistence**: a conversation that outlives the page (servers, headless hosts, a reload). Not a backend to add: the
  storage port went with the durable harness, so it is a redesign; it can start from git
  ([simplification/plan.md](simplification/plan.md) "Later").
