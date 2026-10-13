# Handoff: Prism and the AI packages (2026-10-13)

For the next agent picking up branch `v6-ui-next`. Read this, then [README.md](README.md) (Prism),
[../ai/README.md](../ai/README.md) (AI), and the `AGENTS.md` files (root, `packages/prism`, `packages/ai`). Everything
below is committed and pushed; the working tree was clean.

## Where things stand

**Prism, the new Mol\* UI** (`packages/prism/{ui,viewer,stories,assistant}`): built as mocks in Storybook, with real
runtime pieces underneath. The app wiring (`createPrism`/`createPrismUI`, the plugin behind sources and actions) is not
built yet.

- The skeleton is decided: [design/contexts.md](design/contexts.md) (one mode per app, `PrismContext`/`PrismUIContext`,
  modules, contributions, one command registry with `execute`, origins, history, undo through the mode's checkpointer).
  The spec review ([review.md](review.md)) is resolved up to the "integration step"; the rest of the contracts (review
  M6–M12, M22–M29 details) freeze after it.
- Recent UI work (all in Storybook):
  - **The Inspector** (Shell › Assistant: Trace, Trace running, Trace what it saw, Tools, Tools at a tool, Limits,
    Inspect). It is a `ShellPage`: a full-size card over the canvas, beside the rail and the open panel. Its ways in are
    the Assistant tab's Inspect block, the commands Assistant trace (⌥T), Assistant tools and Assistant limits, and the
    turn card's Trace link.
  - **PD presets in the params form** (Components › Params form › Presets): group presets in the header, flat-group
    presets as a row, object-list presets beside Add, line-graph presets.
  - **PD params as an agent tool's input**: `paramsShape` in `@molstar/prism-assistant/tools/params`.

**The AI packages** (`packages/ai/*`, docs in `.v6/ai`):

- **Ported from Stories Untold**:
  - `agent`, `harness`, `components` (memory and IndexedDB storage, the in-process supervisor; SQLite removed for now);
  - `model-adapters`;
  - `conformance` (the harness suite and `./fibers`);
  - `tools` and `session` ([../ai/tools-and-session.md](../ai/tools-and-session.md));
  - `trace` (feeds the Inspector);
  - `fibers` (durable transcripts);
  - `narration` (TTS providers; audio through Mediabunny with the browser's AAC, then Opus, then WAV; no MP3 encoding,
    because every MP3 encoder is LGPL).
- **Not wired to Prism yet**: no real session feeds the Assistant; everything there is mock data.

## The plan from here

The order was set with the user ([../ai/plan.md](../ai/plan.md) §2;
[../ai/analysis/viewer-mode.md](../ai/analysis/viewer-mode.md) §7):

1. **Next, small**: a **model catalog** in `@molstar/ai-model-adapters` (`./catalog`). Entries are by provider and model
   family, matched by the longest prefix (so a new Gemini version works unchanged and a new family is one entry). Each
   has its label and note, limits, reasoning modes and default, dated prices, and quirks. The adapters' capabilities,
   the session's defaults (replacing `thinkingFor` in `ai-session/src/drive.ts`) and the host's model list and credits
   read it.
2. **Then pause** and write one plan with the Prism UI's ([plan.md](plan.md)), bringing in the sequence viewer (seqstar,
   in `stories-untold-redux/lib/seqstar-next`; Prism draws its own UI over seqstar's headless context:
   [analysis/harness-and-seqstar.md](analysis/harness-and-seqstar.md)). After that, Prism hosts the agent: command
   metadata (`read`, `approval`, `ui`), a selection param with three forms, the core viewer commands, state and focus
   accounts, the checkpointer, `AgentModel` over `AgentSession`, History over fibers, the selections skill.
3. **Decided but not built**:
   - Runtime and lazy contributions (contexts.md §6: `contributions.add/remove`, a stub that loads on first use, with
     video export first). This needs `apps/prism` built as ESM with code splitting ([plan.md](plan.md) Phases 7–8).
   - `ai-toolscript` and `mvs-code` for the story mode.
   - The MCP server.
   - `@molstar/ai-evals`.
   - SQLite storage.
   - Jest as ESM (`.v6/plans/checklist.md`).

Do not start step 2's integration without the user: they want to decide it together.

## How the user works (keep doing this)

- **Style**: flat and quiet, few borders and shadows, the accent only for actions and the picked item, overflowing text
  fades (never "…"), tokens for every color and spacing (`styles:check` enforces the scale), and Mol\*'s own UI ideas
  restyled, not replaced. Check light and dark, wide and narrow, in Storybook, and look before saying it is done.
- **Docs**: short and high-level. Code comments cite design sections (`contexts.md §4`, `harness.md §6`). A change that
  moves a design or plan updates it in the same commit.
- **Commits**: at checkpoints, prefixed `v6 prism:`, `v6 prism-ui:` or `v6 ai:`. The message says what changed and why.
  End it with the `Co-Authored-By` line. Push only when asked; the user often pushes themselves.
- **Files**: new files under `packages/prism` and `packages/ai` carry the Mol\* header plus the line
  `Developed in collaboration with Anyfold <https://anyfold.ai/>` (see any file there).
- **Delegation**: mechanical work (ports, renames, merges) goes to Sonnet subagents. Scope their verification to the
  tests the change touches, never several overlapping full runs: the last SQLite removal took 30 minutes for that
  reason.
- **Opt-in**: prism-ui never imports the modes, the assistant or `@molstar/ai-*`
  (`ui/src/_test/opt-in-boundary.test.ts`). The AI packages never import Mol\*, Prism or rxjs; `node:*` only in
  `ai-components`; zod only in `ai-tools` (Biome enforces these).
- **Decisions the user made, so do not reopen them**:
  - one mode per app, with switching starting a new app;
  - `<prism-viewer>` with a shadow root;
  - the assistant package is `@molstar/prism-assistant`;
  - each mode's agent parts live in an `agent` folder of its package, not separate packages;
  - a denied tool call is told to the model;
  - `ai-components` holds platform code (`src/harness/`, `src/fibers/`);
  - the AI SDK is pinned by major only (everything else is pinned exactly).

## Practical notes

- **Checks**:
  - `npx biome lint .`
  - `npx biome format .`
  - `npx prettier --check "**/*.{md,mdx,yaml,yml,scss}"`
  - `npx tsc -b`
  - `pnpm jest` (excludes `packages/ai`)
  - `pnpm jest:ai`, or `pnpm jest:ai packages/ai/<pkg>` for one package
  - `node scripts/workspace/check.mjs`
  - `pnpm -C packages/prism/<pkg> styles:check`
  - `pnpm -C packages/prism/<pkg> typecheck:storybook`
- **Storybook**: `pnpm -C storybook storybook` (port 6006; the launch config is `storybook`). It lives in the top-level
  `storybook/` package and shows the stories of the packages registered in `storybook/packages.ts`; stories go in
  `_stories/` folders. Another project's Storybook has been holding 6006. For a check, add a temporary launch entry with
  `exec storybook dev -p 6016`, then revert `.claude/launch.json`. Restart Storybook after changing
  `storybook/packages.ts` or `storybook/.storybook/main.ts`.
- **Gotchas**:
  - Every Prism `.tsx` file needs `/** @jsxImportSource @molstar/prism-ui/dom */` first. Without it, Jest compiles it as
    React, and the runtime prints "Cannot render".
  - The Browser pane pauses animations while it is hidden, so measure with timers, not animation frames.
  - The fake-IndexedDB conformance file retries its cases, because fixed-round waits can lose to WebCrypto under load
    ([../ai/plan.md](../ai/plan.md), "The suite waits on conditions").
  - `pnpm-lock.yaml` conflicts are regenerated with `pnpm install`, never hand-merged.
- **Open points to raise with the user**:
  - `thinkingFor` holds Google defaults; the catalog fixes it.
  - The Inspector's mock Trace link always lands on the latest turn.
  - Whether browsers encode AAC and Opus at Gemini's 24 kHz is unverified; there is a runtime check with WAV as the
    fallback.
