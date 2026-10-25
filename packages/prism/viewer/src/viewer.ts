/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { derive } from '@molstar/prism-ui/reactive/state';
import type { AnnotationsApi } from './api/annotations.js';
import { createQueryLanguages } from './targets/languages.js';
import type { ThreadRows } from '@molstar/ai';
import type { PluginContext } from '@molstar/plugin/context';
import type { PluginSpec } from '@molstar/plugin/spec';
import {
  type CanvasInk,
  CommandManager,
  type CommandDefinition,
  createUiTransactionCoordinator,
  type DockId,
  type CommandValues,
  DockIds,
  type DockLayout,
  HistoryManager,
  NavigationManager,
  NotificationsManager,
  observableState,
  type ReadonlyState,
  ViewStateManager,
} from '@molstar/prism-ui/core';
import { type Observable, Subject } from 'rxjs';
import { type ApiModuleKey, type ExtensionPort, type JournalPort, NoJournal, type ViewerApi } from './api/area.js';
import { CoreAreas, createViewerApi, type ViewerApiHost } from './api/index.js';
import type { ViewerCommandContext, ViewerManagers } from './api/managers.js';
import { coreCommands } from './commands/core.js';
import { type PrismViewerOptionValues, resolveOptions } from './options.js';
import { createCheckpointPort, type ViewerCheckpoint, watchOutsideChanges } from './sources/checkpoint.js';
import { mountCanvas, unmountCanvas } from './sources/canvas-host.js';
import { composePluginSpec } from './sources/compose-plugin.js';
import type { ViewerPointerEvent } from './sources/pointer.js';
import { createViewerCamera, type ViewerCamera } from './sources/camera.js';
import { createSequenceSource, type SequenceSource } from './sources/sequence.js';
import { createSceneSequence } from './sources/sequence-scene.js';
import { SequenceApi } from './api/sequence.js';
import type { SequenceBehavior } from './sequence/behaviors.js';
import { DefaultSequenceTracks, type SequenceTrackProvider } from './sequence/tracks.js';
import type { SequenceViewers } from './sequence/model.js';
import type { SavedSequence } from './sequence/runtime.js';
import { SourceScope } from './sources/state.js';
import { createTargetFormat } from './sources/target-format.js';
import { createTargetResolver } from './sources/targets.js';
import { appKeybindings, keymapStore } from './commands/keymap-store.js';
import { isStructureQueryId } from './commands/select-groups.js';
import { watchSelectCommands } from './commands/select.js';
import { loadedStructureData, lociPreviewSink, selectQuerySource } from './sources/select-queries.js';
import { createViewerSources, publicSources, type ViewerSources } from './sources/viewer-sources.js';
import {
  type BuiltInViews,
  isLazyExtension,
  type PrismViewerSequenceSpec,
  type PrismViewerSpec,
  type ViewerExtension,
} from './spec.js';
import { checkSpec, isPartId, PrismViewerSpecError } from './spec-check.js';
import { ExtensionHost, type ExtensionUIPort, isStubDefinition } from './extensions/host.js';
import { Overrides, shellIdOf } from './extensions/overrides.js';
import { createKeptRecordsModel, type KeptRecordsModel } from './journal/kept-records.js';
import {
  createSessionRecord,
  recordFromLines,
  recordLines,
  recordsDatabase,
  type SessionRecord,
} from './journal/record.js';
import { createStorageModel, type StorageModel } from './storage/model.js';
import { viewerStorageSources } from './storage/sources.js';
import { applyExtensionPlugin } from './sources/extension-plugin.js';
import { createLightRestore, type LightRestore } from './sources/restore.js';
import { initDataPalette, paletteLists, setDataPalette } from './sources/palette.js';

// The viewer (viewer-spec.md §6): created from a `PrismViewerSpec`, in a fixed order: validate the spec, compose the
// plugin spec, create and `init()` the plugin, create the sources, managers and API, activate the extensions (lazy ones
// as stubs); a page then mounts the canvas and renders the shell (`createPrismViewer`), headless has neither
// (`createPrismViewerHeadless`). `ready` is set last. `dispose()` undoes it in reverse; a creation that fails undoes
// what it made before throwing.

export type { ViewerPointerEvent } from './sources/pointer.js';

/** The root view models (api.md §11): the shell's on a page; the scene, the selection and the rest as built. */
/**
 * The viewer's view models. The built-in features' (`BuiltInViews`: the Structure panel's `scene`, Home's
 * `homeActions`) are typed as present, as in the default spec; in a viewer whose spec left the feature out, reading one
 * throws, naming the feature (`'scene' in viewer.views` asks).
 */
export interface PrismViewerViews extends BuiltInViews {
  /** The shell (docks, tabs, cards); undefined headless. */
  readonly shell?: import('@molstar/prism-ui').ShellModel;
  /** The one sequence viewer and what it follows (sequence.md); undefined when the spec has none. */
  readonly sequence?: ViewerSequence;
}

/**
 * The sequence view's models: its source (the plugin side, from the start), the canvas ink it is skinned with (set by
 * the page), and the Seq* viewers, made once the spec's runtime is in (sequence.md §2: when a structure is first shown,
 * or an operation or a saved state needs them).
 */
export interface ViewerSequence {
  readonly source: SequenceSource;
  setInk(ink: CanvasInk): void;
  /** The viewers, once loaded; undefined until then. */
  readonly viewers: ReadonlyState<SequenceViewers | undefined>;
  /** Loads the runtime and makes the viewers (once): the viewers and the runtime's module (its views, its readers). */
  load(): Promise<LoadedSequence>;
  /** A behavior (an extension's), attached now or once the viewers are made; returns its removal. */
  addBehavior(behavior: SequenceBehavior): () => void;
  /** The structure tree has each structure's Sequence node (`spec.sequence.tree`, sequence.md §8.1). */
  readonly tree: boolean;
  /** A track provider (an extension's), added now or once the viewers are made; returns its removal. */
  addTrackProvider(provider: SequenceTrackProvider): () => void;
  dispose(): void;
}

export interface LoadedSequence {
  readonly viewers: SequenceViewers;
  readonly runtime: Awaited<ReturnType<PrismViewerSequenceSpec['runtime']>>;
}

export interface ViewerError {
  readonly message: string;
  readonly source: 'plugin' | 'api' | 'extension';
}

export interface PrismViewer {
  /** The engine: an escape hatch, not part of the API's promise (api.md §8). */
  readonly plugin: PluginContext;
  /** Every operation, as a script calls them (origin `script`). */
  readonly api: ViewerApi;
  /** Raw plugin work as one history step (api.md §5); calls it makes through `api` join the step. */
  run<T>(label: string, work: (plugin: PluginContext, api: ViewerApi) => Promise<T>): Promise<T>;
  readonly managers: ViewerManagers;
  readonly views: PrismViewerViews;
  /** Read-only states over the plugin (state.md §4): what views, view models and the agent's account read. */
  readonly sources: ViewerSources;
  /**
   * The camera as a view to follow and set from outside, for keeping viewers' cameras together (`syncCameras`); not
   * the `camera` area's operations, which are a person's or the agent's (api.md §12).
   */
  readonly camera: ViewerCamera;
  readonly state: {
    readonly ready: ReadonlyState<boolean>;
    /** The plugin has been working for a while (its `isBusy`, debounced: a status that does not flicker). */
    readonly busy: ReadonlyState<boolean>;
    /**
     * Something is in progress now: an operation (from its start, a lazy extension's load and a file's download
     * included), an extension loading, or the plugin updating its state. Not debounced: what the empty canvas waits on
     * before it shows the splash, so a slow load never shows it.
     */
    readonly working: ReadonlyState<boolean>;
    /** Nothing loaded. */
    readonly empty: ReadonlyState<boolean>;
  };
  readonly events: {
    readonly hover: Observable<ViewerPointerEvent>;
    readonly click: Observable<ViewerPointerEvent>;
    readonly error: Observable<ViewerError>;
  };
  /**
   * The extensions: which are loaded, each lazy one's status, loading one (README decision 22), and mounting one the
   * spec did not list, at runtime (viewer-spec.md §3).
   */
  readonly extensions: Pick<
    ExtensionHost,
    'loaded' | 'prefixes' | 'status' | 'events' | 'load' | 'ensureFor' | 'mount'
  >;
  /** The session's record (journal.md): what the person and the agent did, in order; in memory until kept. */
  readonly record: SessionRecord;
  /** The records kept in earlier sessions (journal.md §3, §5): Home's list, one opened read-only in History, Forget. */
  readonly keptRecords: KeptRecordsModel;
  /** This tab's light restore, when the spec keeps one (`persist.restore`, journal.md §6): the host runs it. */
  readonly restore?: LightRestore;
  /**
   * What is kept in the browser and Forget (journal.md §7): the viewer's sources (kept records, preferences, the light
   * restore); the assistant registers its own (stored keys).
   */
  readonly storage: StorageModel;
  /** The spec it was made from, and its options resolved. */
  readonly spec: PrismViewerSpec;
  readonly options: PrismViewerOptionValues;
  dispose(): Promise<void>;
}

export interface ViewerCoreOptions {
  readonly spec: PrismViewerSpec;
  /** Makes the plugin from the composed spec (a `HeadlessPluginContext` when images are wanted headless). */
  readonly createPlugin: (spec: PluginSpec) => PluginContext;
  /** What the host offers (api.md §4); `drop` removes one later (no WebGL after all). */
  readonly modules: readonly ApiModuleKey[];
}

/** The viewer without a page: everything but the canvas and the shell. `createPrismViewer` adds those. */
export interface ViewerCore extends PrismViewer {
  readonly host: ViewerApiHost;
  readonly ready$: { set(value: boolean): void };
  /** Undoes what was added since creation (the page's parts), before the core goes. */
  onDispose(fn: () => void | Promise<void>): void;
  /** Registers commands with the spec's keybindings applied (api.md §11: keybindings name commands). */
  registerCommands(defs: readonly CommandDefinition<ViewerCommandContext>[]): () => void;
  /** The spec's overrides, for the shell's own parts. */
  readonly overrides: Overrides;
  /** The page's shell takes the extensions' views. */
  /**
   * Throws for an override of a command nobody registered, and warns for a keybinding of one (viewer-spec.md §4):
   * once every command is in. \`tolerate\`: ids only a page registers (headless).
   */
  checkSpecCommands(tolerate?: (id: string) => boolean): void;
  attachExtensionUI(port: ExtensionUIPort): () => void;
  /** Mounts the canvas in a page's element; false (and the `canvas` module dropped) when it cannot start. */
  mountCanvas(host: HTMLElement): Promise<boolean>;
}

/**
 * What `options.ui.buttons` hides, by option (viewer-spec.md §2.2): the shell's toolbar items by id (the features' by
 * their part ids); `selection-mode` and `trajectory` are the shell's mode switch and trajectory bar.
 */
const ButtonIds: Readonly<Record<string, readonly string[]>> = {
  reset: ['reset'],
  screenshot: ['export'],
  // One button: fullscreen, which v5 calls expand too.
  expand: ['expand'],
  fullscreen: ['expand'],
  settings: ['canvas-settings.toolbar.canvas'],
  selectionMode: ['selection-mode'],
  animation: ['animation.toolbar.animation'],
  trajectory: ['trajectory'],
};

export async function createViewerCore(o: ViewerCoreOptions): Promise<ViewerCore> {
  const { spec } = o;
  checkSpec(spec, Object.keys(CoreAreas));
  const options = resolveOptions(spec.options);
  const eager = (spec.extensions ?? []).filter((e): e is ViewerExtension => !isLazyExtension(e));
  const pluginSpec = composePluginSpec({ base: spec.plugin, extensions: eager, options, finalize: spec.finalize });

  const disposers: (() => void | Promise<void>)[] = [];
  const plugin = o.createPlugin(pluginSpec);
  const undo = async () => {
    for (const d of disposers.splice(0).reverse()) {
      try {
        await d();
      } catch (e) {
        console.error(e);
      }
    }
    plugin.dispose();
  };

  try {
    await plugin.init();

    const made = createViewerSources(plugin);
    disposers.push(() => made.dispose());
    const sources = made.sources;

    const empty = sources.empty;
    const commandContext = observableState<ViewerCommandContext>({ empty: empty.value });
    const emptySub = empty.subscribe((e) => commandContext.set({ empty: e }));
    disposers.push(() => emptySub.unsubscribe());

    // This viewer's own query languages: the page-wide defaults, then the spec's (targets.md §1; README decision 29).
    const languages = createQueryLanguages(spec.languages ?? []);
    const targetFormat = createTargetFormat(plugin, languages);
    const history = new HistoryManager<ViewerCheckpoint>(createCheckpointPort(plugin), {
      checkpoints: options.ui.undoSteps,
    });
    const notifications = new NotificationsManager();
    const { viewState, navigation } = layoutManagers(spec, options);
    const transactions = createUiTransactionCoordinator();
    const managers: ViewerManagers = {
      commands: new CommandManager<ViewerCommandContext>({
        context: commandContext.asReadonly,
        formats: { get: (name) => (name === 'target' ? targetFormat : undefined) },
        // The command menu's highlighted rows light their atoms (the Select group's queries, commands.md §5).
        preview: lociPreviewSink(plugin),
        // The spec's keys over the commands' own; the person's remaps, kept with the preferences (commands.md §6).
        keymap: { app: appKeybindings(spec), store: keymapStore(spec.persist?.preferences) },
      }),
      transactions,
      history,
      navigation,
      viewState,
      notifications,
    };
    disposers.push(() => {
      history.dispose();
      transactions.dispose();
      notifications.retire();
    });
    disposers.push(watchOutsideChanges(plugin, history));

    const targets = createTargetResolver(plugin, languages);
    const modules = new Set(o.modules);
    // Operations that load what an extension owns ask the extensions first; set once they are added.
    const extensionPort: ExtensionPort = { ensureFor: async () => {} };
    // The record's part of session files (journal.md §5): bound once the record exists, below.
    const journalPort: JournalPort = { ...NoJournal };
    const host = createViewerApi({
      plugin,
      sources,
      targets,
      languages,
      modules: { has: (k) => modules.has(k) },
      managers,
      extensions: extensionPort,
      journal: journalPort,
    });
    disposers.push(() => host.dispose());

    const keys = new Map((spec.ui?.keybindings ?? []).map((k) => [k.command, k.key]));
    // Every command goes through here: the core's, the shell's, extensions' and stubs' (overrides.ts).
    const overrides = new Overrides(spec.overrides);
    // The spec's keybindings are the keymap's app layer (above), not rewritten into the definitions.
    const registerCommands = (defs: readonly CommandDefinition<ViewerCommandContext>[]) =>
      managers.commands.register(overrides.commands(defs));

    // The extensions, in order (viewer-spec.md §3): eager ones activated, lazy ones as stubs.
    const person = host.api.as({ kind: 'person' });
    const viewScope = new SourceScope();
    disposers.push(() => viewScope.dispose());
    const sequenceSpec = spec.sequence || undefined;
    const sequence = sequenceSpec
      ? createSequence(viewScope, plugin, host, o.modules, sequenceSpec, overrides, (e) =>
          error.next({
            message: `The sequence view did not load: ${e instanceof Error ? e.message : e}`,
            source: 'extension',
          }),
        )
      : undefined;
    if (sequence) disposers.push(() => sequence.dispose());
    disposers.push(registerCommands(coreCommands(person)));
    // The Select group: a command per selection query, made when wanted after the structures change (commands/select.ts).
    disposers.push(
      watchSelectCommands({
        queries: selectQuerySource(plugin),
        api: person,
        targets,
        structures: sources.structures,
        data: () => loadedStructureData(plugin),
        registerLazy: (source) =>
          managers.commands.registerLazy({
            make: () => {
              const made = source.make();
              return made && { commands: overrides.commands(made.commands), groups: made.groups };
            },
          }),
        registerGroups: (groups) => managers.commands.registerGroups(groups),
      }),
    );
    const extensions = new ExtensionHost({
      host,
      person,
      registerCommands,
      openCommand: (id) => void managers.commands.open(id),
      annotations: () => host.api.area('annotations') as AnnotationsApi | undefined,
      ...(sequence
        ? { addSequenceBehavior: (b) => sequence.addBehavior(b), addSequenceTrack: (t) => sequence.addTrackProvider(t) }
        : {}),
      applyPlugin: (parts) => applyExtensionPlugin(plugin, parts),
      overrides,
      ...(spec.persist?.preferences ? { preferences: spec.persist.preferences } : {}),
      post: (text, kind) => notifications.post(text, kind),
      task: (label) => {
        const id = `extension:${label}`;
        notifications.putTask({ id, label });
        return () => notifications.removeTask(id);
      },
    });
    disposers.push(() => extensions.dispose());
    // One record per viewer session (journal.md §1), in memory; History reads it, the person may keep it, in the spec's
    // rows when it gives them (`persist.rows`, viewer-spec.md §2.6; typed loosely there, a `ThreadRows` here).
    const keepIn = spec.persist?.rows as ThreadRows | undefined;
    const record = await createSessionRecord({
      api: host.api,
      undone: history.undone,
      redone: history.redone,
      ...(keepIn ? { keepIn } : {}),
    });
    disposers.push(() => record.dispose());
    // The kept records after a reload (journal.md §3, §5): Home lists them, History opens one read-only.
    const keptRecords = createKeptRecordsModel({
      record,
      ...(keepIn ? { rows: keepIn } : {}),
      database: recordsDatabase(spec.persist?.key),
      max: options.ui.keptRecordsMax,
    });
    disposers.push(() => keptRecords.dispose());
    journalPort.lines = () => recordLines(record);
    journalPort.show = (name, lines) => {
      const file = recordFromLines(lines);
      if (!file) return false;
      keptRecords.show({ id: `file:${name}`, title: name, first: file.first, last: file.last }, file.reader);
      return true;
    };
    // What is kept and where, with Forget (journal.md §7): the viewer's own sources; the Storage view lists them.
    // The origin that holds it, where there is one (a file: page's is 'null').
    const origin = globalThis.location?.origin;
    const storage = createStorageModel(origin && origin !== 'null' ? { origin } : {});
    const own = viewerStorageSources({
      ...(spec.persist ? { persist: spec.persist } : {}),
      record,
      resetKeybindings: () => managers.commands.keymap.resetAll(),
    });
    for (const s of own.sources) disposers.push(storage.register(s));
    disposers.push(() => {
      own.dispose();
      storage.dispose();
    });
    // The data palette the options ask for (accessibility.md): before anything is loaded, so no step. Its state is
    // registered either way, so a session's palette comes back.
    initDataPalette(plugin, paletteLists(spec.palettes));
    const palette = options.representation.defaultPalette;
    if (palette !== 'classic') await setDataPalette(plugin, palette);
    // The light restore (journal.md §6), written after each step when the spec gives it a store.
    const restore = spec.persist?.restore
      ? createLightRestore(viewScope, plugin, history, spec.persist.restore, {
          notify: (message) => notifications.post(message, 'info'),
          // late-bound: the port is wired once the extensions are added (below)
          ensureFor: (need) => extensionPort.ensureFor(need),
        })
      : undefined;
    // Features behind an option (viewer-spec.md §3.1: remote states) are left out unless it is on.
    const gated: Readonly<Record<string, boolean>> = { 'remote-states': options.ui.showRemoteStates };
    const listed = (spec.extensions ?? []).filter((e) => gated[e.id] !== false);
    await extensions.add(listed);
    extensionPort.ensureFor = (need) => extensions.ensureFor(need);
    // Each structure's sequence in the structure tree (sequence.md §8.1), through the API as the person.
    const sceneSequence =
      sequence && sequenceSpec?.tree !== false
        ? createSceneSequence(
            viewScope,
            sequence.source,
            sequence.viewers,
            person,
            () => void managers.commands.open('sequence.legend'),
          )
        : undefined;
    // The built-in features' view models (viewer-spec.md §3.1), made once their areas are there; headless too.
    const builtInViews: Partial<BuiltInViews> = {};
    for (const e of listed) {
      if (isLazyExtension(e) || !e.builtIn || !e.views) continue;
      Object.assign(
        builtInViews,
        e.views({
          scope: viewScope,
          plugin,
          person,
          openCommand: (id, values) => void managers.commands.open(id, undefined, values as CommandValues),
          ...(sceneSequence ? { sequence: sceneSequence } : {}),
        }),
      );
    }
    // An override of a command nobody registers is a mistake (viewer-spec.md §4); view parts are checked by spec-check.
    // An override of a command nobody registers is a mistake (viewer-spec.md §4), and a keybinding of one is said:
    // checked once every command is in, the page's after its shell (`checkSpecCommands`), headless at once.
    const checkSpecCommands = (tolerate: (id: string) => boolean = () => false) => {
      // The Select families of the loaded structures exist only once one is loaded: a well-formed id is known.
      const unknown = overrides
        .unused()
        .filter((k) => !isPartId(k) && !tolerate(k) && !isStructureQueryId(k) && !managers.commands.get(k));
      if (unknown.length) {
        throw new PrismViewerSpecError(unknown.map((k) => `override '${k}': there is no command '${k}'`));
      }
      for (const id of keys.keys()) {
        if (!tolerate(id) && !isStructureQueryId(id) && !managers.commands.get(id)) {
          console.warn(`ui.keybindings: there is no command '${id}' (its key is not bound).`);
        }
      }
    };

    // What the person sees in the status column: results of the person's and the URL's runs, failures, the plugin's
    // errors, running tasks, undo turning off for a large scene.
    const error = new Subject<ViewerError>();
    const subs = [
      // An extension that did not load, or whose views did not (the page also says it in the status column).
      extensions.events.failed.subscribe((f) =>
        error.next({ message: `${f.id} could not load: ${f.message}`, source: 'extension' }),
      ),
      host.runner.finished.subscribe((r) => {
        const result = r.result;
        if (!result || r.nested || r.group) return;
        const def = host.runner.definitions.value.get(r.operation);
        if (!result.ok) {
          if (result.code === 'aborted' || result.code === 'confirm') return;
          if (r.origin.kind === 'person' || r.origin.kind === 'url' || r.origin.kind === 'link') {
            notifications.post(result.hint ? `${result.message} ${result.hint}` : result.message, 'warning', {
              detail: durationOf(r.ms),
            });
          }
          error.next({ message: result.message, source: 'api' });
        } else if (
          (r.origin.kind === 'person' || r.origin.kind === 'url') &&
          def &&
          !def.read &&
          // View operations (hover's highlight, focus) say nothing; a lazy stub's run is the real one's result.
          (def.history !== false || isStubDefinition(def))
        ) {
          // A quiet one (a story's step: the card shows it) is said to screen readers only.
          if (def.quiet) notifications.announce(result.did);
          // How long it took: in the Log only, never on the canvas (journal.md §3).
          else notifications.post(result.did, 'info', { detail: durationOf(r.ms) });
        } else if (r.origin.kind === 'person' && def && !def.read && /^(select|focus|measure)\./.test(r.operation)) {
          // A selection, a focus or a measurement by target is how a person without sight picks (accessibility.md): its
          // result is said (the live region), with no message on show.
          notifications.announce(result.did);
        }
      }),
      plugin.events.log.subscribe((e) => {
        if (e.type === 'error') {
          notifications.post(String(e.message), 'error');
          error.next({ message: String(e.message), source: 'plugin' });
          return;
        }
        // The state's own lines, one per transform it created or updated, with its time (v5's log): the Log's
        // Transforms, within the Log's one bound. Never on the canvas.
        const transform = TransformLine.exec(String(e.message));
        if (transform) {
          // `12ms`, `1.204s` (core's `formatTimespan`) as the runs' are said: `12 ms`, `1.204 s`.
          const detail = transform[2].replace(/^([\d.]+)(ms|s)$/, '$1 $2');
          notifications.log(transform[1], 'info', 'Transforms', { detail });
        } else if (e.type === 'warning') {
          notifications.log(String(e.message), 'warning', 'Plugin');
        }
      }),
      sources.tasks.subscribe((tasks) =>
        notifications.setTasks(
          tasks.map((t) => ({
            id: String(t.id),
            label: t.label,
            ...(t.progress !== undefined ? { progress: t.progress } : {}),
            ...(t.overlay ? { overlay: true } : {}),
            ...(t.abort ? { abort: t.abort } : {}),
          })),
        ),
      ),
      history.checkpointingOff.subscribe(() =>
        notifications.post('This scene is large: undo is off. Settings turns it back on.', 'warning'),
      ),
    ];
    disposers.push(() => {
      for (const s of subs) s.unsubscribe();
      error.complete();
    });
    notifications.mount();
    disposers.push(() => notifications.unmount());

    const ready = observableState(false);
    const viewerCamera = createViewerCamera(plugin);
    disposers.push(() => viewerCamera.dispose());
    let disposed = false;
    const core: ViewerCore = {
      plugin,
      api: host.api,
      run: (label, work) =>
        host.runner.opaque(label, { kind: 'script' }, (context) =>
          work(plugin, host.bind({ origin: { kind: 'script' }, within: context })),
        ),
      managers,
      views: viewsOf({ ...builtInViews, ...(sequence ? { sequence } : {}) }),
      // Read-only (api.md §6): the marks' writers are the API's (`highlight.mark`).
      sources: publicSources(sources),
      camera: viewerCamera.camera,
      state: {
        ready: ready.asReadonly,
        busy: sources.busy,
        working: derive(
          [sources.updating, host.runner.running, managers.notifications.tasks],
          (updating, running, tasks) => updating || running.length > 0 || tasks.length > 0,
        ),
        empty,
      },
      events: { hover: sources.pointer.hover, click: sources.pointer.click, error: error.asObservable() },
      extensions,
      record,
      keptRecords,
      ...(restore ? { restore } : {}),
      storage,
      spec,
      options,
      host,
      ready$: ready,
      registerCommands,
      overrides,
      checkSpecCommands,
      attachExtensionUI: (port) => extensions.attachUI(port),
      mountCanvas: async (el) => {
        const ok = await mountCanvas(plugin, el, { pauseOutOfView: options.ui.pauseOutOfView });
        // Disposed while it mounted: the disposers have run, so the canvas is let go here.
        if (disposed) {
          if (ok) unmountCanvas(plugin);
          return false;
        }
        if (ok) disposers.push(() => unmountCanvas(plugin));
        else modules.delete('canvas');
        return ok;
      },
      onDispose: (fn) => disposers.push(fn),
      async dispose() {
        if (disposed) return;
        disposed = true;
        ready.set(false);
        await undo();
      },
    };
    return core;
  } catch (e) {
    await undo();
    throw e;
  }
}

/** The layout's and navigation's starting state from the spec (viewer-spec.md §2.5, §2.2). */
/** A state transform's log line (`Created Model in 12ms.`): what it did, and how long it took. */
const TransformLine = /^((?:Created|Updated) .+) in (\S+)\.$/;

/** A run's duration as the Log says it: "120 ms", "2.4 s". */
function durationOf(ms: number | undefined): string | undefined {
  if (ms === undefined) return undefined;
  return ms < 1000 ? `${Math.round(ms)} ms` : `${(ms / 1000).toFixed(ms < 10_000 ? 1 : 0)} s`;
}

function layoutManagers(spec: PrismViewerSpec, options: PrismViewerOptionValues) {
  const layout = spec.ui?.layout ?? {};
  const docks: Partial<Record<DockId, Partial<DockLayout>>> = {};
  const nav: Partial<Record<DockId, { active?: string; collapsed?: boolean }>> = {};
  for (const id of DockIds) {
    const d = layout.docks?.[id];
    if (!d) continue;
    const { tabs, size, floating, pinned, active, collapsed } = d;
    // Tabs by part id (`core.tab.home`, `snapshots.tab.snapshots`) or the core's short ids.
    docks[id] = {
      ...(tabs ? { tabs: tabs.map(shellIdOf) } : {}),
      ...(size !== undefined ? { size } : {}),
      ...(floating !== undefined ? { floating } : {}),
      ...(pinned !== undefined ? { pinned } : {}),
    };
    nav[id] = {
      ...(active !== undefined ? { active: shellIdOf(active) } : {}),
      ...(collapsed !== undefined ? { collapsed } : {}),
    };
  }
  // A tab placed by the spec: listed in its dock, by order, so it lands there when it registers.
  const placed = Object.entries(layout.tabs ?? {}).sort((a, b) => (a[1].order ?? 0) - (b[1].order ?? 0));
  for (const [id, place] of placed) {
    const tab = shellIdOf(id);
    const dock = place.dock ?? 'right';
    const current = docks[dock]?.tabs ?? [];
    docks[dock] = { ...docks[dock], tabs: [...current.filter((t) => t !== tab), tab] };
  }
  const buttons: Record<string, boolean> = {};
  for (const [option, on] of Object.entries(options.ui.buttons)) {
    if (!on) for (const id of ButtonIds[option] ?? [option]) buttons[id] = false;
  }
  const viewState = new ViewStateManager({
    docks,
    ...(layout.mode ? { mode: layout.mode } : {}),
    compactRail: !!layout.compactRail,
    buttons,
  });
  const navigation = new NavigationManager({ docks: nav });
  return { viewState, navigation };
}

/**
 * The sequence (sequence.md): its source and its area in the API from the start, its viewer over the structure
 * the source follows once the runtime is in, and its state in snapshots (`prism-viewer.sequence`, sequence.md §7.1),
 * kept as given until then.
 */
function createSequence(
  scope: SourceScope,
  plugin: PluginContext,
  host: ViewerApiHost,
  modules: readonly ApiModuleKey[],
  spec: PrismViewerSequenceSpec,
  overrides: Overrides,
  failed: (e: unknown) => void,
): ViewerSequence {
  const link = host.api.as({ kind: 'link' });
  const source = createSequenceSource(scope, plugin, link);
  const ink = observableState<CanvasInk>({ theme: 'light', color: '#ffffff' });
  const viewers = observableState<SequenceViewers | undefined>(undefined);
  const initial = spec.style;
  // Before the runtime: the saved state given (a snapshot, a session), applied once it is in; an extension's behaviors.
  let pending: { value: unknown } | undefined;
  const behaviors = new Map<SequenceBehavior, (() => void) | undefined>();
  const tracks = new Map<SequenceTrackProvider, (() => void) | undefined>();
  let disposed = false;
  let loading: Promise<LoadedSequence> | undefined;
  let loaded: LoadedSequence | undefined;
  const load = (): Promise<LoadedSequence> =>
    (loading ??= (async () => {
      const runtime = await spec.runtime();
      if (disposed) throw new Error('The viewer is gone.');
      const made = runtime.createSequenceViewers({
        ink: ink.asReadonly,
        structure: source.structure,
        actions: source.actions,
        dataByUrl: source.dataByUrl,
        scheduler: runtime.sequenceFrames(modules.includes('page')),
        // The viewer's own behaviors and track providers, overrides applied; the extensions' come as they activate
        // (sequence.md §3).
        behaviors: overrides.parts(spec.behaviors),
        tracks: overrides.parts(spec.tracks ?? DefaultSequenceTracks),
      });
      const viewer = made.get();
      if (initial === 'text') viewer.setStyle('compact');
      for (const b of behaviors.keys()) behaviors.set(b, viewer.addBehavior(b));
      for (const t of tracks.keys()) tracks.set(t, viewer.addTrackProvider(t));
      if (pending) runtime.applySaved(viewer, (ref) => source.setFollowed(ref), pending.value, initial);
      pending = undefined;
      // A loaded SeqViewSpec (sequence.md §6, David 2026-10-15) shows while the sequence shows its structure: set, the
      // source shows that structure; once the person focuses or loads another (the source ends `show`), it ends.
      scope.add(viewer.override.subscribe((o) => source.show(o?.structure)));
      scope.add(
        source.shows.subscribe((ref) => {
          if (ref === undefined && viewer.override.value) void viewer.setOverride(undefined);
        }),
      );
      loaded = { viewers: made, runtime };
      viewers.set(made);
      return loaded;
    })().catch((e) => {
      loading = undefined;
      throw e;
    }));
  // A structure to show: the runtime loads (its chunk beside the data's download).
  scope.add(
    source.structure.subscribe((s) => {
      if (s && !loading) load().catch(failed);
    }),
  );
  scope.add(host.addArea('sequence', (ctx) => new SequenceApi(ctx, load, source)));
  scope.add(
    plugin.state.registerCustomState('prism-viewer.sequence', {
      get: (): SavedSequence | undefined => {
        // A structure that is gone leaves its edits out (sequence.md §7.1).
        const exists = (ref: string) => plugin.state.data.cells.has(ref);
        if (loaded) return loaded.runtime.savedOf(loaded.viewers.get(), source.followed.value, exists);
        // Not loaded: what was given, as given (nothing given: the defaults, which a saved `undefined` means).
        return pending?.value as SavedSequence | undefined;
      },
      set: (value) => {
        if (loaded) loaded.runtime.applySaved(loaded.viewers.get(), (ref) => source.setFollowed(ref), value, initial);
        else {
          pending = { value };
          // What it follows is the source's: now.
          const v = value as SavedSequence | undefined;
          source.setFollowed(v && 'structure' in v ? v.structure : undefined);
        }
      },
    }),
  );
  return {
    source,
    tree: spec.tree !== false,
    setInk: (i) => ink.set(i),
    viewers: viewers.asReadonly,
    load,
    addBehavior(b) {
      const made = viewers.value;
      behaviors.set(b, made ? made.get().addBehavior(b) : undefined);
      return () => {
        behaviors.get(b)?.();
        behaviors.delete(b);
      };
    },
    addTrackProvider(t) {
      const made = viewers.value;
      tracks.set(t, made ? made.get().addTrackProvider(t) : undefined);
      return () => {
        tracks.get(t)?.();
        tracks.delete(t);
      };
    },
    dispose() {
      disposed = true;
      viewers.value?.dispose();
    },
  };
}

/** Each built-in view model's feature, for the error that names it (`PrismViewerViews`). */
const BuiltInViewFeatures: Readonly<Record<keyof BuiltInViews, string>> = {
  scene: 'StructurePanelFeature',
  presets: 'StructurePanelFeature',
  look: 'StructurePanelFeature',
  quickStyles: 'StructurePanelFeature',
  homeActions: 'HomeFeature',
};

/** The views object: a built-in one the spec left out throws when read, naming its feature. */
function viewsOf(present: Partial<PrismViewerViews>): PrismViewerViews {
  return new Proxy(present as PrismViewerViews, {
    get(target, prop, receiver) {
      if (typeof prop === 'string' && !(prop in target) && prop in BuiltInViewFeatures) {
        const feature = BuiltInViewFeatures[prop as keyof BuiltInViews];
        throw new Error(`The '${prop}' view model is not in this viewer: add ${feature} to its spec's extensions.`);
      }
      return Reflect.get(target, prop, receiver);
    },
  });
}
