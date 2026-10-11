# The harness and the sequence viewer: what Prism has to host

Read-only sweeps (2026-10-11) of anystar (the agent harness, `stories-untold/lib/anystar/*`; the same code in
`stories-untold-redux`) and seqstar (the sequence viewer, `stories-untold-redux/lib/seqstar` and the rewrite in
`lib/seqstar-next`, spec in `designs/seqstar/spec`). What each is, how Stories Untold hosts it, and what that asks of
Prism's contracts before they are frozen ([../review.md](../review.md)).

## 1. anystar: a turn runner, not a command registry

- **Packages** (all UI-free, no zod): `@anystar/agent` (messages, `ToolDefinition` with a JSON Schema `inputSchema`, the
  `ModelAdapter` port), `@anystar/harness` (the host protocol, `createRunner`, the drive, tools, nested runs,
  observation), `@anystar/harness-components` (storage: memory, IndexedDB, sqlite; an in-process supervisor),
  `@anystar/model-adapters` (Anthropic, OpenAI responses/chat/compatible, Google, over the `ai` SDK),
  `@anystar/harness-conformance` (tests).
- **The loop.** Prompts are submitted as items with an idempotency key; a runner holds a lease and commits each step (a
  model call or a batch of tool calls) atomically to durable storage, so a turn survives a reload and resumes. Policy is
  a "drive" of decision points (permission, schedule, budget, failure, model-result, repair, boundary, interrupt,
  next-turn). Streaming is `Observation` events. Stop is durable (`client.stop`); each tool call gets an `AbortSignal`.
- **What it needs from a host:** storage, `ModelBindings` (keys stay in the host: sealed in the browser, env on a
  server), `ToolBindings` (`resolve(definition) → invoke(ctx): AsyncIterable<ToolFrame>`), a clock and ids, a drive, and
  `configure(items)` (model, tools, limits).
- **How Stories hosts it** (`lib/stories/stories-agent`, `stories-application`): tools are defined once (zod, turned
  into JSON Schema, arguments loosely normalized first), marked `read` (parallel, de-duplicated) or write (operations on
  the document, all or nothing), `approval` for risky ones, and the capability `modules` they may use. Every prompt
  carries a compact state block and the person's focus. A whole turn is one undo step, attributed to the agent; user
  edits are locked while a turn runs; history is an inverse/redo list, capped and persisted. The command menu's
  assistant is a one-off read-only harness that proposes a filled command. An MCP server (`apps/server/src/mcp.ts`,
  streamable HTTP, loopback) offers the same tools; a batch of calls is one edit.

**What Prism's contracts need** (review M5, M22–M25):

1. **Commands as tools.** A command is the tool: stable id, description, and its params as JSON Schema (a
   `ParamDefinition` → JSON Schema converter, plus a loose argument normalizer). Metadata: `read`, `approval`,
   capabilities (and so whether an agent outside the app may call it), and grouping for deferred loading.
2. **Programmatic execution.** `execute(id, args, { origin, signal })` returning a result the harness can carry: status,
   a model-readable summary, `transient`, and `effect: none | applied | unknown` (so retries are safe); progress while
   it runs (the tool frames).
3. **Origin on every change.** `person | agent | mcp | behavior | link`, carried from the command to the history and to
   the sources' events (see §2: seqstar tags changes the same way to skip echoes).
4. **History and undo.** One step per person action and per agent turn; one checkpoint kept until the next action
   (decided); while a turn runs, the person's edits wait or run outside the turn (decide per mode: Stories locks them).
5. **State and focus accounts.** `PrismContext` gives the prompt a state block and the focus (what is selected, picked,
   shown), per mode.
6. **The assistant port is async** (`propose(request, signal)`), a one-off read-only run.
7. **Headless.** The harness, the tools and an MCP server run without a page: commands and their contexts must run in
   `PrismContext` alone (review M2), and the agent package keeps the harness, the tool layer and the adapters out of the
   UI.

## 2. seqstar: a headless viewer with its own spec

- **What it is.** A reusable sequence viewer: tracks and text layouts, many layer types (letters, ranges, points,
  variants, alignments with trees, logos, matrices, arcs), numbering schemes, streaming to genome scale, SVG/PNG export.
- **SeqViewSpec (SVS 1.0)**: a JSON document shaped like MolViewSpec (`{kind, ref, params, custom, children}`, refs
  mandatory), data kinds (sequence, ligand, complex, alignment, tree, mapping, annotation, …) and view kinds (view,
  section, track, layer, focus). Ref-addressed patches with inverses; snapshots are titled patches over one document. A
  structure link is optional (`custom.mvs` on units).
- **Architecture (seqstar-next).** One headless `SeqViewerContext` (`createSeqViewer(spec, ports)`), states as rxjs
  observables with a synchronous `.value`, Canvas2D drawing through a render backend, DOM only for chrome. Hosts mount
  it with `mountSeqViewer(container, ctx, ui)` or `<seq-viewer>`, give it ports (loaders, a document port that decides
  what is committed), and talk to it through states, events and JSON commands (`seqv:select`, `seqv:focus`, …). Theming
  via `--seq-*` tokens; `chrome: "clean" | "full"`.
- **Coupling to Mol\*: none in code.** The designed link is bidirectional (hover, selection, focus), over `Coordinates`
  and registered coordinate systems for the structure's chains; changes carry an `origin` to skip echoes.
  `@seqstar/molstar` is specified, not built.
- **Status.** The rewrite's models, spec, renderers and behaviors are done; its canvas surface and DOM/React UIs are not
  (only the old `lib/seqstar/viewer` mounts today).

**What Prism's contracts need** (review B1, M18):

1. **A mode attaches a model that owns an outside context**: the viewer mode's sequence model creates the
   `SeqViewerContext` with Prism's loaders and document port, and disposes it with the app (the module slot of B1).
2. **A view that mounts a foreign renderer** (`elementView` or a mount function) into a dock tab, sized by its host.
3. **The structure link is Prism's to write**: registering a structure's chains as seqstar systems, translating loci to
   `Coordinates` and back, over the selection, highlight and focus sources, with origins to skip echoes.
4. **Stories**: a story document can carry SVS next to MVS (both ref-addressed, patchable, snapshot-able); the story
   model treats both as parts of one document.
5. **Theming**: map `--msp-*` onto `--seq-*` in the tab's scope; `chrome: "clean"` for a story's reader, `"full"` for
   the viewer's dock.
6. **Until the rewrite mounts**, the build keeps a placeholder sequence tab (decided); the contracts above are what it
   plugs into.
