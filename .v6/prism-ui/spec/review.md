# Spec review (2026-10-15)

A multi-agent review of the spec as committed in `159bccf09`: six lenses (consistency, claims versus code, architecture,
v5 parity, agent and security, implementability) found 80 issues, merged into 54; two skeptics per finding kept 45 as
confirmed (both agreed) and 4 as plausible (one agreed), and refuted 5. Then a report and a completeness critic. To be
folded into the spec; each fix cites its finding id.

## How the findings were folded in

Fixed the same day, with David's answers to the report (README decisions 24–32). Where a finding led to a decision
rather than the edit it proposed, the decision is named.

| Findings                         | Where                                                                                                   |
| -------------------------------- | ------------------------------------------------------------------------------------------------------- |
| F1, F4, F5, F10, F34             | sequence.md §7.1–7.3 (checkpoints through `plugin.state`, `HistoryPort`, groups, runs in flight, failures, changes outside the API, `registerCustomState`); operations.md `history?: false`; decision 24 |
| F23, F35                         | sequence.md §3 (stable refs, removals with their dependents, a composition test, the person's edits overlay), §7.1 (`{ version, style, colors, structure, rail, edits }`) |
| F2, F11, F12                     | agent.md §1, §3 (`createViewerAgent`, `paramsShape` in the viewer, a dynamic toolkit); viewer-spec.md §2.4, §3 (`tools`); decision 25 |
| F3, F20                          | README §3.5 (framing by camera target, not `viewOffset`; canvas sizing after the WebGPU work); decision 28 |
| F6, F44, F45                     | viewer-spec.md §3 (stubs declare every part and operation), §3.1 (one module per feature); README §3.4 (lazy config, format order, conflicts) |
| F7, F8, F25, F38, F40            | one id scheme (viewer-spec.md §4), operation ids aligned (sequence.md §5), History everywhere, statuses and references, managers in library.md |
| F9, F30, F33                     | parity.md §6, §6a; decision 27 (extensions register commands; model-archive and ANVIL designed anew; annotations as a command) |
| F13, F14, F15, F16, F17          | agent.md §4 (per-address approval with domain rules, undo-dependent preauthorization, quoted data, key table, the origin); app.md §1; decision 26 |
| F18, F19, F52                    | journal.md §1 (an index of kept records), §2 (heartbeat while parked), §6 (when the light restore runs) |
| F21                              | sequence.md §4–5 (the link's writes go through operations; the mapping's direction)                     |
| F24                              | [api.md](api.md) §4 (`ApiContext`, `ApiModuleKey`, `Projection`, `ViewerApi`)                          |
| F26, F27                         | app.md §4 (the one URL list), §2 (`FullPrismViewerSpec()`, an artifact of the app); viewer-spec.md §7 |
| F28, F31, F32, F49, F50, F51     | state.md, parity.md, operations.md (additions history, keyboard, selection layers and Modify, `load.snapshot`, `highlight` and hover/click events, focus entries) |
| F29                              | operations.md `load` (core `files`, `action`; the loaders feature the rest)                             |
| F36, F37, F42, F43, F46, F48     | viewer-spec.md §2.5, README §3.1, the `pluginExtension` example, README §4, operations.md §1, library.md §3 |
| F39                              | targets.md §2                                                                                           |
| Critic 1                         | decision 29, library.md §3, viewer-spec.md §2.6; viewports and several state trees in README §6 (step 2) |
| Critic 2                         | operations.md §1 (identity, all-or-nothing failure), agent.md §3 (`ask`); decision 30                   |
| Critic 3                         | README §3.7; decision 31                                                                                |
| Critic 4                         | journal.md §7; decision 32                                                                              |
| Critic 5                         | journal.md §2, §6 (summaries past a size, skipped restore writes); budgets per scale stay to be measured |

## 1. History, checkpoints and undo

**F1 · major · confirmed** The checkpoint the spec describes ("the plugin's snapshot through the snapshot manager") has side effects, captures the wrong state, and crashes when there is no canvas.
- Evidence:
  - `snapshots.ts` `getStateSnapshot` → `syncCurrent` adds or replaces entries, and `setStateSnapshot` → `clear()`.
  - `state.ts:335-346` defaults are `behavior: false`, `structureSelection: false`, `camera: true`. Restore resets focus and animates the camera.
  - `state.ts:89-99` dereferences `canvas3d!.camera`.
  - Selection mode and snapshot entries are not in `Snapshot`.
  - viewer-spec §6 activates loaders, which run URL calls, before the canvas mounts.
- Edit, sequence.md §7:
  - Define the checkpoint port as `plugin.state.getSnapshot({ data: true, behavior: false, camera: false, animation: false, image: false, structureSelection: <decide>, … })` and restore with `plugin.state.setSnapshot`. Never go through `managers.snapshot`.
  - List what undo restores (selection mode, the sequence's choices and snapshot entries are captured beside the snapshot) and what it never touches (camera, hover, the Snapshots list; strip `structureFocus`).
  - State what undo across a lazy extension's load does.
  - Define the 5 MB size measure (exclude RawData), or budget by time only.
- Edit elsewhere:
  - Fix README §3.3 in both places ("through its snapshot manager"; "what a plugin snapshot holds").
  - viewer-spec §6: URL calls run after the canvas mounts or after `ready`.
  - Add a headless test: undoable op, then Undo, with no canvas.

**F4 · major · confirmed** History semantics are undefined in six places:
- how a turn is grouped into one entry;
- what happens to person steps made during a turn (Undo vs Undo turn, and which entry is "latest");
- Undo while a run is in flight;
- redo (listed, but entries hold only a before-checkpoint; there is no `canRedo` and no journal kind);
- which ops are steps (`select`/`focus`/`camera`);
- state.md §3 still says snapshots are taken only "before an agent's turn".

Evidence: README §3.3 says "Undo restores the latest". journal.md §3 has "Undo turn". operations.md §2 lists `history: undo, redo`. README §3.6 `managers.history` has no `canRedo`. contexts.md §5 had a release rule that the spec dropped.
- Edit, sequence.md §7:
  - Type `HistoryPort { checkpoint(); restore(c) }` and give the manager `beginGroup(origin)` / `endGroup()`, keyed by `origin.turn`.
  - Interleaving rule: a person's call closes the turn's group, or it waits for the turn.
  - Undo turn is enabled only while the turn is the latest entry, or it confirms before discarding later steps.
- Edit, operations.md:
  - Undo is disabled during a run, or aborts it through `signal` and then restores.
  - Drop redo, or define it fully (redo stack, `canRedo`, a journal kind).
  - State that `select.*`, `focus.*`, `camera.*` and `view.*` are not steps.
- Fix state.md §3.

**F5 · major/minor · confirmed** "A change outside the API ends any held agent checkpoint" has no detector.
- Behaviors write the tree constantly: the focus representation on every click, volume streaming on camera moves.
- The spec itself writes outside the API: Seq* load cells (sequence.md §2) and SVS follow into 3D (§6). A tree-diff restore deletes those cells.
- A restore is itself an outside change.

Evidence: `structure-focus-representation.ts:192-218,265`; `volume-streaming/behavior.ts:384`; `core/state/state.ts:169-172` (setSnapshot = updateTree).
- Edit, sequence.md §7, add "what counts as inside":
  - runs under `api.run`/`viewer.run`;
  - the history's own restore (suppressed);
  - behavior-owned cells, allow-listed by transformer or tagged.
  - Any other change invalidates every entry it makes unsafe, not only the agent checkpoint.
- Edit, sequence.md §2/§6: route Seq* loads and SVS loads through operations (origin `link`, plus `load.svs`), or tag their cells so a restore keeps them.
- Update README §3.2 to point here.

**F10 · major/minor · confirmed** `Snapshot.custom` as a plain field gives the viewer no hook. The manager's add, syncCurrent and session save call `getSnapshot()`, and apply, applyNext and play call `setSnapshot()` internally, so `'prism-viewer.sequence'` is never written or applied.
- Evidence: `state.ts:87-122` builds a fixed object literal; `snapshots.ts:166,255,273`.
- Edit, sequence.md §7 and README §4:
  - Add `PluginState.registerCustomSnapshotState(key, { get(); set(v | undefined) })`.
  - `getSnapshot` fills every registered key. `setSnapshot` calls `set` for every key after the tree is applied, passing `undefined` (reset) when the key is absent.
  - Unknown keys pass through unchanged.
  - Transition frames do not apply `custom`.

**F23 · major/minor · confirmed** "Behaviors patch disjoint parts, no merge rule" fails under SeqViewSpec's all-or-nothing, no-cascade patches.
- Other behaviors' tracks reference the complex that follow owns, so follow's remove-and-insert on a structure switch fails with `svs.ref.unresolved`.
- Person edits ("patches too") have no layer. §7 saves only `{style, theme, structure, rail}`, so hidden tracks are lost on undo, snapshot apply and session open.

Evidence: `sequences/_spec/format.md:426-440`; `viewer.md:72-87`; sequence.md:72-76, 128-136.
- Edit, sequence.md §3:
  - Follow keeps stable refs and switches with `set`, never remove-and-insert.
  - A behavior that removes a node also removes the nodes that refer to it, in the same patch.
  - Add a cross-behavior composition test.
- Edit, person edits: define a viewer-owned overlay keyed by track ref, applied after the behaviors.
- Edit, sequence.md §7: add `edits` to `'prism-viewer.sequence'` and to checkpoints.

**F34 · minor · confirmed** "View tools change only the view (focus, select)" contradicts "Three layers of state", which puts selection and focus in engine state. operations.md has no flag for "changes state, not a step", and selection mode is not in `Snapshot`.
- Evidence: README §3.3; `state.ts:392-393`; `context.ts:133`.
- Edit:
  - operations.md §1: add `history?: false` and apply it to `view.*`.
  - Reword View tools (README §3.3, agent.md §2) to "never the structures' content".
  - Say selection mode sits outside the plugin snapshot.

## 2. Agent toolkit and session composition

**F2 · major · confirmed** Extension skills cannot reach the toolkit.
- `ToolGroup.tools` holds ids only. Registrations go separately in `ToolkitSpec.tools`. `createToolkit` validates once and throws on an unregistered id.
- No spec field carries an extension's `ToolRegistration`s.
- The lazy stub (`{id, title, summary}`) contradicts "every tool is declared from the start".
- Registry enums "refreshed" after creation hit the same wall.
- A stub has no operation ids for `api.run`.
- (Mitigation: `session.ts` `catchUp` can append definitions between turns.)

Evidence: `ai/tools/src/toolkit.ts:45-63,80-89,143-180`; `session.ts:217,332-335`.
- Edit, viewer-spec §3 and §2.4: add `tools?: readonly ToolRegistration<ViewerApi>[]`.
- Edit, lazy extensions, pick one and record it in README §3.3 and agent.md §3:
  - (a) the stub declares full definitions with schemas generated at build time and proxy handlers that load the chunk; or
  - (b) loading rebuilds the session between turns via `catchUp`, and the spec states what that costs the cache.
- Edit, operations.md §2: registry args are open strings checked at call time, or a refresh is a rebuild.
- Edit, stub: add `operations?`, so `api.run` loads the chunk and replays the call.
- Note that `list_skills` is a viewer tool.

**F11 · major/minor · confirmed** The spec has two owners for the AgentSession. journal.md §4 and agent.md §1 say the lazy, DOM-bound prism-assistant creates it. README §4 says the viewer composes ai-session. Headless evals and MCP need the same session and policies, and `AssistantFactory` is undefined.
- Evidence: README §4 vs journal.md §4; viewer-spec §2.4; `prism/assistant/package.json`.
- Edit:
  - Put a DOM-free `createViewerAgentSession(viewer, { model, policies })` in the viewer's `agent/`, exported from main and `/headless`.
  - The assistant only renders it.
  - Define `AssistantFactory = (viewer, session) => contribution`, with the assistant depending on the viewer and never the reverse.
  - agent.md §6: MCP and evals use the same factory.
  - journal.md §3: HistoryTab moves to the viewer.

**F24 · major/minor · confirmed** Core types are undefined or named twice:
- `ApiContext` vs `ViewerApiContext`;
- `ApiModule` and `ApiModuleKey`;
- what `ViewerApi` is: a toolkit module record (keyed by area? by host capability?) that also carries `run`, `running`, `started`, `finished`;
- `Projection`, which collides with ai-tools' exported `Projection`;
- `ViewerSources`;
- the per-slot context for `ViewRenderer<Ctx>`.

Evidence: operations.md:11-22; viewer-spec.md:139-150; `ai/tools/src/registry.ts:25`; `prism/ui/src/dom/view.ts:23`.
- Edit: add a "Core types" block to operations.md §1, referenced from viewer-spec §3 and sequence.md §3:
  - one `ApiContext` and what it holds;
  - `ApiModule`;
  - `ViewerApi` with an explicit key set, used the same way in README §3.3 and `Operation.modules`;
  - rename `Projection` to the prism-ui optimistic contribution type;
  - one `Origin` union;
  - a slot context map.

**F12 · minor · confirmed** `paramsShape` lives only in `@molstar/prism-assistant`, which the viewer may not depend on. No spec file says it moves, and the AI specs still cite the old path. Its dependencies are all reachable from the viewer, so this is a move, not a blocker.
- Evidence: `prism/assistant/src/tools/params.ts:42`; `ai/_spec/tools.md:58`; `ai/_spec/README.md:90`; the AI plan (since removed).
- Edit:
  - agent.md §1: `paramsShape` and its test move to `packages/prism/viewer/src/agent/params.ts`.
  - Update the three AI-spec references.

## 3. Agent security

**F13 · minor/major · confirmed** The `approval` tag is static per operation, but the network policy depends on the arguments ("URLs the person wrote"). The spec never maps this to the toolkit's `approvedBy(args, api)` and never defines "the person wrote", so URLs from tool results or metadata are not excluded. Documents the agent loads (MVS, sessions, SVS) fetch every URL they contain. Latent cases: `settings.set` on endpoint servers, and `state.update` on Download cells.
- Evidence: agent.md §4; operations.md §1; `config.ts:57-63`; `ai/_spec/tools.md:37`.
- Edit, agent.md §4:
  - Preauthorized = changes only viewer state and fetches only known-service hosts.
  - Network approval uses `approvedBy`: https, plus a known service or a URL verbatim in the person's own message. Never URLs from results or documents.
  - A loaded document is checked by every URL it names.
- Edit, operations.md: mark endpoint settings and data-source state edits `ui-only` or `approval: 'network'`.

**F15 · major/minor · confirmed** The untrusted-text rule covers tool results only. File-controlled labels (titles, MVS labels and tooltips, snapshot names, cell labels) go into the per-turn state and focus blocks, which sit inside the user message next to the request.
- Evidence: agent.md §4, §5; `ai/_spec/session.md` Composition; `ai/session/src/compose.ts`.
- Edit:
  - agent.md §4/§5: data-derived labels go into state and focus as quoted, length-capped string fields (for example JSON), and the guidance says they are data.
  - Add an eval fixture: an MVS label that contains an instruction.

**F14 · minor · confirmed** Preauthorization rests on "the history undoes the turn". A large-scene turn over the hard limit runs with no checkpoint, and the turn card only informs. `export.*` and `settings.set` fall in neither class, and `undoable: false` is not tied to `approval`.
- Evidence: agent.md §4; sequence.md §7; operations.md §1.
- Edit:
  - agent.md §4: writes in a turn with no checkpoint ask, or the card asks once before the turn runs.
  - Place `export.*` and `settings.set` explicitly.
  - `undoable: false` requires an `approval`.
  - Give the hard limit a number.

**F16 · minor · confirmed** Row 1 of the key table claims a non-extractable CryptoKey in IndexedDB defeats "copying storage". It does not: a profile copy holds both the key and the ciphertext.
- Evidence: agent.md §4 table.
- Edit:
  - Change "Against" to "reading one store alone; not a profile copy".
  - Show "scrambled, not protected" in every row.

**F17 · minor · confirmed** Stores are origin-scoped, and molstar.org is shared by every Mol* app, demo and the docs. Any of those pages can use the key (non-extractable still allows use) and read kept records. The spec only excludes "code running in the page".
- Evidence: agent.md §4; app.md §1; `.github/workflows/docs.yml`.
- Edit:
  - agent.md §4 and journal.md §5: "or on its origin".
  - app.md §1: no third-party or untrusted pages on that origin, or give the app its own subdomain.
  - viewer-spec §2.6: an embedder's `persist` stores under the host page's origin.

## 4. Extensions, lazy loading, packaging

**F6 · major · confirmed** Several mechanisms need part metadata that lives in unloaded code:
- The lazy stub type covers only commands, sections, skills and `owns`. README promises "a toolbar item".
- `viewer.state-tree` (a tab) cannot get a rail button before load, and `viewer.remote-states` cannot show the setting that enables it.
- §6 validates overrides and duplicates before lazy parts are known.
- Lazy plugin conflicts arrive after the "throw at init(), atomically" point.

Evidence: viewer-spec §3 stub type; §3.1 default list; §4, §6; README §3.4.
- Edit:
  - The stub declares every part id the extension will register (tabs with title, icon and dock; toolbar; rows; settings; API modules), and a load that registers an undeclared id fails.
  - Validate at creation against eager parts plus stub declarations.
  - An eager extension's `ui()` resolves at mount.
  - README §3.4: lazy conflicts fail that extension's load atomically and are reported through `failed`.

**F29 · major · confirmed** Loading is fixed core in viewer-spec §3.1 (Home, drag and drop, `.molx`), but operations.md gives the whole `load` area to the optional `viewer.loaders`. The only other path for registry actions is the lazy `state.action`.
- Evidence: viewer-spec §3.1; operations.md §2 `load` row; parity.md §1.
- Edit: split the area.
  - `load.files` and `load.action` (a registry root StateAction by id) are core.
  - The loaders feature adds `entry`, `url`, `mvs`, `snapshot` and the URL parameters.
  - Update README §3.2, parity rows 15/17/62, and state.md rows 40-41.

**F9 · major/minor · confirmed** Moving model-archive out of prism-viewer conflicts with the viewer's "own" pLDDT and PAE behaviors and with `load.entry afdb`. Keeping them gives a package cycle; moving them leaves pLDDT colouring and the AFDB preset broken in the default spec. Separately, rcsb and pdbe are "data-only" but add sequence tracks.
- Evidence: `prism/viewer/src/sources/viewer-plugin.ts:12,141,203-215,266`; `extensions/plugin/src/loaders.ts:160-179`; `scripts/workspace/check.mjs:233`.
- Edit, pick one:
  - (A) pLDDT/PAE behaviors, the theme, predicted colouring and the AFDB preset move to model-archive's `/prism-viewer` entry; or
  - (B) the viewer refers to them by id and hides them when the provider is absent (§3.4 Slim viewers).
- Fix sequence.md §3 and parity.md §6 to match.
- Give rcsb and pdbe a `/prism-viewer` entry, or let `pluginExtension` take `sequence`.

**F27 · minor · confirmed** The single-file example uses `DefaultPrismViewerSpec()`, which holds built-ins only. The file exports nothing that refers to the extension list, so the bundler drops it. (Packaging is fine: the distribution already builds private apps by path.)
- Evidence: app.md §2; viewer-spec §5; `scripts/distribution.mjs:73-74`.
- Edit:
  - app.md §2: name the entry `apps/prism-viewer/src/single-file.ts`, exporting `MolstarPrismViewerSpec()` (built-ins plus extensions).
  - Use it in the example.
  - Cross-reference it from viewer-spec §3.1.

**F44 · minor · confirmed** "Exported both ways" (`ViewerFeatures.x` and `ViewerFeatures.lazy.x`) on one runtime object pulls the eager code into the entry chunk. Stub location is unstated: in the extension's entry it loads the heavy code, written in the app it drifts from `owns`.
- Evidence: viewer-spec §3.1; README §3.1; app.md §1.
- Edit:
  - Lazy forms live in modules that never import the eager ones.
  - Each lazy extension ships a light `/prism-viewer-stub` entry.
  - Add tests: stub vs loaded extension, and a bundle check that the entry chunk holds no lazy modules.

**F45 · minor · confirmed** Three gaps in the lazy plugin parts:
- Their `config` is unspecified, and items read in the constructor (e.g. `HistoryCapacity`) would be ignored.
- Lazy formats register after the custom formats, so `auto()` ties resolve differently.
- Duplicate behaviors throw one at a time after `register()`, which is not atomic.

Evidence: `context.ts:96-97,618-620`; `transient.ts:326`; `formats/registry.ts:67-77`.
- Edit, README §3.4:
  - No `config` in lazy plugin parts (narrow the type).
  - Lazy formats need an explicit `priority`.
  - `compose-plugin.ts` checks duplicate behavior ids before `init()` and lists them all.

**F40 · minor · confirmed** README §3.4 says the plugin part is overridable per extension, but viewer-spec §4 has no plugin part kind and no id for one. (The commands/transactions half of the original finding is not a real issue.)
- Edit: in viewer-spec §4, add `<ext>.plugin` as a part kind, or say plugin parts are changed by passing `{ ...ext, plugin }`.

**F42 · minor · confirmed** The example passes a raw transformer where `PluginSpec.Behavior` is required.
- Evidence: `spec.ts:94-104`; `apps/viewer/src/extensions.ts:56`.
- Edit, viewer-spec §3 and README §3.4: `behaviors: [PluginSpec.Behavior(RCSBValidationReport)]`.

**F30 · minor · confirmed** anvil is listed as data-only, but its "Membrane Orientation" generic representation control needs a row. state.md says `genericRepresentationControls` is unused.
- Evidence: `extensions/anvil/src/behavior.ts:56`.
- Edit, pick one:
  - give anvil a small `/prism-viewer` entry; or
  - define a viewer-owned generic row over `genericRepresentations` and update state.md.

**F43 · minor · plausible** Importing `@molstar/plugin-extension/loaders` adds plugin-ui and React peers to prism-viewer's dependencies. This costs install size only; the bundle is unaffected.
- Edit, README §4: record this as an accepted cost, or split the loaders into a UI-free package.

**F46 · minor · plausible** "/spec has no plugin import" really means type-only. Built-in features need plugin value imports outside `api/` and `sources/`, and nothing enforces the boundary for prism-viewer.
- Edit:
  - viewer-spec §1: "`import type` only".
  - Name the folder for feature plugin parts.
  - Add the import check that contexts.md §7 promises.

## 5. Canvas and styling

**F3 · major · confirmed** `camera.viewOffset` is already the jitter channel: `setJitter` overwrites it every sample and `clearJitter` disables it, and the renderer feeds it ×16 into `uViewOffset`. `paddedFocus` never uses viewOffset, so the cited precedent does not exist. shell.md and ShellView shrink the canvas for the pinned panel and the rail, and that change is not in §4.
- Evidence: `canvas3d/passes/jitter.ts:102-112`; `renderer.ts:611-616`; `prism/viewer/src/sources/framing.ts:13-27`; `shell.module.css:51`.
- Edit, README §3.5 step 1:
  - Add a persistent projection inset to Canvas3D's camera, composed under jitter and honoured by picking, screenshots and the passes.
  - Drop the `paddedFocus` analogy.
  - "Steps 1–4 are Canvas3D changes."
- Edit, §4: add graphics and prism-ui rows (ShellView stops insetting the canvas).
- Update shell.md §2.

**F20 · major/minor · confirmed** The plugin already calls `requestResize`, and Canvas3D never sizes the canvas element itself. The debounced cost is `resizeCanvas`/`setCanvasSize`. Setting `canvas.width` in a ResizeObserver callback clears the buffer, and the draw waits for the next frame, giving one blank frame per drag step.
- Evidence: `context.ts:379-383,405-412`; `canvas3d.ts:814-818,1411`; `canvas3d/util.ts`.
- Edit, README §3.5 step 2:
  - The resize event only records a pending size. Element sizing and `syncPixelScale` move into Canvas3D's `render` before `passes.updateSize`.
  - Update the §4 row.
  - "Stretched" should read "kept at its old size".

**F48 · minor · confirmed** `adoptStyles` takes only a Document and MspElement passes `ownerDocument`, so a shadow-root `<prism-viewer>` renders unstyled. Token defaults on `:root` match nothing inside the shadow root. styling.md still says "no per-shadow-root injection".
- Evidence: `prism/ui/src/styles/registry.ts`; `dom/element.ts:65`; `styles/theme.css:15`.
- Edit, library.md:
  - `adoptStyles(Document | ShadowRoot)`, with every adopted root updated on `defineStyles`.
  - MspElement uses `getRootNode()`.
  - Add `:host` next to `:root`.
  - Overlays mount inside the viewer's root.
- Add a test with two elements and a lazy extension loaded after mount. Mark styling.md:18 as superseded.

## 6. Journal and persistence

**F18 · major · confirmed** Four problems with the light restore:
- It and the loaders' URL actions both run at start, with no order or winner, so reloading `?pdb=1tqn` loads twice.
- sessionStorage survives same-tab navigation, so editing the query restores the old scene.
- "From the byte cache" cannot happen: that cache is in-memory AssetManager state and is empty after a reload.
- Download cells with headers or body (credentials) are written to storage.

Evidence: journal.md §6; app.md §4; viewer-spec §5/§6; `core/src/util/assets.ts:20-27,90`.
- Edit, journal.md §6:
  - Restore only when the navigation type is `reload` and the stored `location.search` matches; in that case skip the URL actions, otherwise drop the restore.
  - Say "HTTP cache".
  - Leave out cells with headers or body.
- Edit, viewer-spec §6: put both steps in the creation order (after mount, before `ready`).

**F19 · minor · confirmed** Kept records cannot be listed (RowStore and FiberHost have no listing), and "copies into IndexedDB and keeps appending" has no mechanism.
- Evidence: `ai/fibers/src/storage.ts:29-51`; `host.ts:84`.
- Edit, journal.md §1/§5:
  - The viewer keeps its own catalog (id, title, times), updated on keep and forget, which Home and the bound of 20 read.
  - Define the copy: replay settled entries, re-begin an open turn, switch the writer.
  - The bound is advisory across tabs.
  - README §4 note: the catalog is the viewer's.

**F52 · minor · confirmed** A turn parked for approval keeps an idle fiber stream. Inactivity (hard bound 1 h) aborts it, and the later `settle` fails.
- Evidence: `ai/_spec/threads.md` (then `fibers.md`) inactivity and heartbeat; `ai/fibers/src/host.ts:573`.
- Edit, journal.md §2:
  - Settle at park with outcome `awaiting-approval` and begin a continuation entry on resume, or heartbeat open turns.
  - State `inactivityTimeoutMs`.

**F25 · minor · confirmed** journal.md §3 makes History replace Log. viewer-spec §3.1/§2.2/§7, library.md, parity.md §4 and state.md §1 still have a Log tab and `showLog`.
- Edit:
  - Replace Log with History (part id `viewer.tab.history`) in all of them.
  - `showLog` becomes `showHistory` or is dropped, and `layoutShowLog` maps accordingly.

## 7. Sequence and targets

**F21 · minor · confirmed** `coordinatesOf` goes 3D→Seq* (`link.fromLoci` runs → Seq* Coordinates). Seq*→3D is `link.toLoci`/`toMembers`, and no loci→author-address step exists yet. `linkViews` defaults to `focus: true`, a camera move outside the API that the spec never mentions.
- Evidence: `sequences/molstar-poc/src/viewer.ts:524-536,640,671`; `plugin.ts:124`.
- Edit:
  - targets.md §3 and sequence.md §5: `link.toLoci`, then targets' projection to author addresses.
  - Note the camera focus as accepted engine interaction, or pass `focus: false`.

**F35 · minor · confirmed** The sequence's state fields are named four ways:

| Source | Fields |
|---|---|
| README §3.6 | `style/theme/structure/legend` |
| sequence.md §7 | `style/theme/structure/rail` |
| operations.md | `style/colors/follow/legend` |
| state.md | "style, colors, structure followed" |

- Edit:
  - Define the record once in sequence.md §7 (`style`, `colors` or `theme`, `structure`, `rail`, and `legend` if saved).
  - Use those names everywhere.
  - Add a `rail` op.

**F39 · minor · confirmed** `substructureFromRows` is placed "in the MVS runtime", but it is in `@molstar/model` schema-selection. Following the citation leads to plugin-dependent `@molstar/mvs`. (The any-molql module-name half of the original finding is wrong.)
- Edit, targets.md §2: cite `StructureElement.Schema.toExpression` and `schema-selection/row-structure.ts`.

## 8. Naming and consistency

**F7 · minor · confirmed** sequence.md §5 and parity.md use ops and areas that operations.md §2 does not define:
- `view.show`, `structure.color`, `selection.set`, `focus(target)`, `navigation.reveal` (a manager, bypassing the API);
- `structure.contents.*`, `selection.*`, `preset.*`, `viewState`, `zenodo.*`.
- Edit:
  - Rewrite them to `focus.set`/`camera.frame`, `representation.color`, `view.open`, `select.*`, `contents.*`, `structure.preset`, `view.*`, plus the extension's fixed id.
  - operations.md §4: the parity column names §2 areas only.

**F8 · minor · confirmed** Part ids use four schemes:
- viewer-spec §4: `<owner>.<kind>.<name>`;
- README: `geo-export.export-section`, `kinemage.views`;
- sequence.md: `sequence.follow-focus` vs §4's `viewer.sequence.follow`;
- built-in extension ids such as `viewer.export.image`, which parse as owner/kind/name.
- Edit:
  - One rule in §4, with `viewer` as a reserved owner and dot-free built-in extension ids.
  - Fix the examples in README §3.1, sequence.md §3 and operations.md.

**F26 · minor · confirmed** The spec states three ways who reads URL parameters (README: the loaders feature, renderer options included; viewer-spec and app.md: the options PD). Some parameters have no target field:
- `pick-padding` (not even in v5 config);
- `map-provider` → `volumeStreamingServer` (unstated);
- `debug-mode`/`timing-mode` (app-level);
- `show-toggle-fullscreen`.

Also missing:
- URL parameters: `structure-url-format`, `-is-binary`, `pdb-dev`;
- v5 options: `layoutShowLeftPanel`, `collapseRightPanel`, `rcsbAssemblySymmetry*`.

parity.md cites `triggerLabel`, but the field is `applyLabel`.
- Edit:
  - README §3.2: loaders read load parameters only; the options PD reads the rest.
  - Add the mappings and the missing entries; assembly-symmetry options go to its extension's settings.
  - Rename to `applyLabel`.

**F36 · minor · confirmed** `ui.layout: ShellOptions & …` brings in `ShellOptions.buttons`, a second way to hide buttons alongside "one way to hide".
- Evidence: `prism/ui/src/shell/shell.ts:111-127,431`.
- Edit: `Omit<ShellOptions, 'buttons'>`, and say the viewer drives `setButtons` only from overrides.

**F37 · minor · confirmed** README says "an extension cannot add a rail button", but `ViewerExtensionUI.tabs` adds exactly that. README's slot list also differs from `ViewerSlot`.
- Edit:
  - Reword: tabs come only through `ui.tabs`, placed by the viewer and `ui.layout.tabs`.
  - Align the slot list.

**F38 · minor · confirmed** Housekeeping:
- viewer-spec §2.2 cites a nonexistent `api.md`.
- journal.md has §6 before §5.
- parity.md is "draft" while README says the neighbouring files are accepted.
- Edit:
  - Point to app.md §4.
  - Reorder journal.md.
  - Align the status lines, or note parity.md as a live draft.

## 9. Parity gaps

**F33 · minor · confirmed** v5's default "Automatic (w/ Annotation)" preset (pLDDT, QMEAN, SB-NCBR) has no owner, so `?afdb=` opens with plain colouring.
- Evidence: `apps/viewer/src/presets.ts`; `plugin-spec.ts:107`; `app.ts:70`.
- Edit:
  - Add a parity row.
  - Name the owner and set `DefaultRepresentationPreset` in config.
  - Resolve provider ids at runtime with an AutoPreset fallback, so there is no model-archive dependency.

**F28 · minor · confirmed** state.md calls `additionsHistory` "undo last selection". It is the ordered last-24 picks that Measure (top 2/3/4) and Superpose › atoms read; the toolbar Undo is `state.data.undo()`.
- Evidence: `selection.ts:41,235`; `measurements.tsx`; `selection.tsx:256`.
- Edit:
  - Fix state.md §1 and §3.
  - operations.md: measure and superpose take ordered target lists that default to the picks.
  - Map toolbar Undo to `history.undo`.

**F49 · minor · confirmed** `load.snapshot` is named in README §3.2 and app.md §4 runs `snapshot-url`/`snapshot-id`, but the `load` row lacks it.
- Edit:
  - Add `load.snapshot` (url plus type, or remote id), working whatever the remote-states setting.
  - Mark the other loaders script-only.

**F50 · minor · confirmed** parity.md has no row for `Viewer.structureInteractivity` or for hover/click subscriptions. Embedders can still use `viewer.plugin`.
- Edit: add a parity row mapping them to `applyStructureInteractivity(viewer.plugin, …)` and `viewer.plugin.behaviors.interaction.*`.

**F51 · minor · confirmed** The focus picker (per-ligand, branched-residue and BIRD entries split by operator) and focus history are not covered.
- Evidence: `plugin/ui/src/structure/focus.tsx`.
- Edit: split the parity row, and give each part a place or mark it later or dropped.

**F31 · minor · plausible** The parity row "keyboard (none beyond Enter/Escape)" is wrong for canvas and behavior keys. Coexistence is already settled in design/commands.md §6, and keys reach only the viewer under the pointer.
- Edit: reword the row to list the trackball, behavior and snapshot keys, "kept; conflicts shown in Help (commands.md §6)".

**F32 · minor · plausible** The selection theme layers (emissive, material, clipping, resets) and component Modify-by-Selection have no ops. The lists are marked first pass.
- Edit:
  - Add `representation.emissive`/`material`/`clipping` with reset (or one `layer` op) and `component.modify`.
  - Fix parity.md's `selection.*`.

## Edits to make first

1. sequence.md §7: write the checkpoint port, the history type, the interleaving/redo/run rules, what counts as "inside", and the custom-snapshot registry (F1, F4, F5, F10, F23).
2. viewer-spec §3 and agent.md: add `tools` to extensions, decide how lazy skills declare tools, define `AssistantFactory`, and put the DOM-free session factory in the viewer (F2, F11, F24, F12).
3. README §3.5: replace `viewOffset` with a camera inset and move canvas sizing into `render` (F3, F20).
4. Lazy stubs: declare every part, validate against the declarations, and say where stubs live (F6, F44, F45).
5. operations.md `load`: split core and feature, add `load.snapshot` and `load.action` (F29, F49).
6. journal.md §6: add the restore gate and creation-order slot, fix the byte-cache claim (F18).
7. agent.md §4: per-call network approval through `approvedBy`, quoting of data-derived labels, fixes to the key-table claims and an origin note (F13, F15, F16, F17).
8. Settle model-archive and the default annotation preset (F9, F33).
9. One pass to make the names agree: operation ids, part ids, sequence state fields, Log→History (F7, F8, F25, F35).

## Completeness critic

1. **Several viewers on one page, and viewers inside a host page.** Look at README §2.12 ("several on a page"), app.md §3, viewer-spec.md §2.5–2.6 and app.md §4. The spec never says which viewer gets a key press. Today's code listens on the whole page. `packages/prism/ui/src/testing/shell.tsx:368` adds a `window` keydown listener for ⌘K/⇧⌘P and every command shortcut. `views/shell/canvas-chrome.tsx:209` adds one on `document`. So two `<prism-viewer>`s both react to every key, and an embedded viewer takes the host page's ⌘K and ⌘Z. Other shared state:
   - `persist.key` is a single namespace (`'prism-viewer'`). Two instances with `persist` would overwrite each other's sessionStorage restore and write into the same IndexedDB rows.
   - app.md §4 treats `debug-mode` and `timing-mode` as per-viewer options, but they map to `setDebugMode` and `setTimingMode` in `packages/core/src/util/debug.ts:44,48`, which are process-wide globals.
   - Browsers cap the number of live WebGL contexts. The spec has no plan for a page that opens many viewers, and doesn't mention context loss, although Canvas3D already has `contextLost`/`contextRestored`.

   No finding covered key scope, store namespacing or global setters.

2. **What a run leaves behind when it fails or is aborted, and what "desired state" means for operations that create things.** Look at operations.md §1 ("Never throws", "Desired state", the single run path with `signal`) and README §3.3. A step opens with a checkpoint taken before it. The spec doesn't say whether a failed or aborted run restores that checkpoint, or leaves a half-built state tree. Loaders commit several times (download, parse, model, structure, preset), so a failure partway leaves partial state. Idempotence is defined only for setters (`representation.set`). For `load.entry`, `component.add`, `measure.distance`, `snapshot.add` and `superpose.*`, "twice changes nothing" needs a rule for recognising that the thing already exists. One example: does a second "load 1TQN" do nothing, or load a second copy? A retried agent step is only safe if this is defined. F4 covers history ordering; nothing covered partial failure or identity for operations that create things.

3. **Decision 15 (mobile, touch, screen reader) has no place in any spec file.** The only text is README §2.15. Nothing describes:
   - A path to canvas picking, measuring and focusing that doesn't need sight, such as the scene tree plus typed targets.
   - Whether `ApiResult.did` and the agent's streamed reply are announced. Today there is a single `aria-live` (`canvas-chrome.tsx:458`).
   - What replaces hover-only features on touch: the hover label, sequence hover, hover cards.
   - How measurements and superposition, which need ordered picks, work with touch.
   - Reduced motion for README §3.5's pinning, which is animated in the WebGL render loop. styling.md's reduced-motion reset only finishes CSS animations.

   No lens checked accessibility, yet the decision says these are required.

4. **Versioning of everything stored or deployed.** Look at app.md §1, journal.md §6, sequence.md §7 and viewer-spec.md §8.2.
   - `'prism-viewer.sequence'` in `Snapshot.custom` has no version.
   - Preferences and keybindings in localStorage are keyed by part and command ids, which §8.2 says will be named as the code is written.
   - Because `index.html` is not cached, a reload after a deploy reads the old sessionStorage light restore with new code (transformer params may have changed).
   - app.md says a version's files are immutable. A GitHub Pages publish (`peaceiris/actions-gh-pages` in `.github/workflows/docs.yml`) normally replaces the branch contents. A tab left open across a deploy would then fail on its next lazy `import()` of an old hashed chunk. That's the extension `failed` status, with no recovery described.

   No finding covered migration of stored data or chunk lifetime across deploys.

5. **Cost and scale outside the renderer.** sequence.md §7 sets a budget only for history checkpoints, yet the same serialization runs elsewhere with no budget or off switch:
   - The light restore serializes the snapshot "after each step (debounced)" (journal.md §6).
   - Every resolved target returns "residues by author number per chain" (targets.md §2), and that list goes into results, the agent's account and each `run` record entry. For `all` on a ribosome or a mesoscale scene, that is very large.
   - Scene and contents records, the per-turn state block, and the `dataByUrl` index rebuilt from cell events (sequence.md §2) have no stated limits for scenes with thousands of components.

   No performance lens was used apart from the canvas resize (F20).

## Refuted

- F22: history: false on the document port does not turn off Seq*'s undo
- F41: focusBehavior changes the base behavior list, not params set at the options layer
- F47: Public extension packages would depend on Prism packages that are private today, and the spec does not say how
- F53: Command links in file-supplied markdown run commands with person origin, so no approval applies
- F54: No first buildable slice is named, though the spec's dependencies force an order
