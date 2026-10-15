/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { FiberRows } from '@molstar/ai-fibers';
import type { PluginContext } from '@molstar/plugin/context';
import type { PluginSpec } from '@molstar/plugin/spec';
import {
  type CanvasInk,
  CommandManager,
  type CommandDefinition,
  createUiTransactionCoordinator,
  type DockId,
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
import type { ApiModuleKey, ExtensionPort, ViewerApi } from './api/area.js';
import { CoreAreas, createViewerApi, type ViewerApiHost } from './api/index.js';
import type { ViewerCommandContext, ViewerManagers } from './api/managers.js';
import { coreCommands } from './commands/core.js';
import { type PrismViewerOptionValues, resolveOptions } from './options.js';
import { createCheckpointPort, type ViewerCheckpoint, watchOutsideChanges } from './sources/checkpoint.js';
import { mountCanvas, unmountCanvas } from './sources/canvas-host.js';
import { composePluginSpec } from './sources/compose-plugin.js';
import type { ViewerPointerEvent } from './sources/pointer.js';
import { createHomeActions } from './sources/home.js';
import { createSceneSource } from './sources/scene.js';
import { createSequenceSource, type SequenceSource } from './sources/sequence.js';
import { SequenceApi } from './api/sequence.js';
import { DefaultSequenceTheme } from './sequence/document.js';
import { createSequenceViewers, type SequenceStyle, type SequenceViewers } from './sequence/model.js';
import type { TracksRail } from './sequence/presets.js';
import { rafScheduler, type Scheduler } from '@molstar/seq-viewer';
import { SourceScope } from './sources/state.js';
import { createComponentOptionsForm, createPresetsModel, createQuickStylesModel } from './sources/structure-tab.js';
import { createTargetFormat } from './sources/target-format.js';
import { createTargetResolver } from './sources/targets.js';
import { createViewerSources, type ViewerSources } from './sources/viewer-sources.js';
import { isLazyExtension, type PrismViewerSpec, type ViewerExtension } from './spec.js';
import { checkSpec, isPartId, PrismViewerSpecError } from './spec-check.js';
import { ExtensionHost, type ExtensionUIPort, isStubDefinition } from './extensions/host.js';
import { Overrides, ShellCommandIds, shellIdOf } from './extensions/overrides.js';
import { createSessionRecord, type SessionRecord } from './journal/record.js';
import { createStorageModel, type StorageModel } from './storage/model.js';
import { viewerStorageSources } from './storage/sources.js';
import { applyExtensionPlugin } from './sources/extension-plugin.js';
import { createLightRestore, type LightRestore } from './sources/restore.js';
import { setDataPalette } from './sources/palette.js';
import type { HomeAction } from './model/home.js';
import type { PresetsModel, QuickStylesModel, SceneModel } from './model/scene.js';
import type { ParamsFormModel } from '@molstar/prism-ui/forms/params-form';

// The viewer (viewer-spec.md §6): created from a `PrismViewerSpec`, in a fixed order: validate the spec, compose the
// plugin spec, create and `init()` the plugin, create the sources, managers and API, activate the extensions (lazy ones
// as stubs); a page then mounts the canvas and renders the shell (`createPrismViewer`), headless has neither
// (`createPrismViewerHeadless`). `ready` is set last. `dispose()` undoes it in reverse; a creation that fails undoes
// what it made before throwing.

export type { ViewerPointerEvent } from './sources/pointer.js';

/** The root view models (README §3.3): the shell's on a page; the scene, the selection and the rest as they are built. */
export interface PrismViewerViews {
  /** The shell (docks, tabs, cards); undefined headless. */
  readonly shell?: import('@molstar/prism-ui').ShellModel;
  /** The scene tree and inspector (the Structure tab). */
  readonly scene: SceneModel;
  /** The representation presets, from the registry (the Structure tab's Preset menu). */
  readonly presets: PresetsModel;
  /** The component manager's options as a form (the Structure tab's options). */
  readonly componentOptions: ParamsFormModel<any>;
  /** Quick styles (the Structure tab's Styles). */
  readonly quickStyles: QuickStylesModel;
  /** Home's loading actions, from the registry. */
  readonly homeActions: ReadonlyState<readonly HomeAction[]>;
  /** The one sequence viewer and what it follows (sequence.md); undefined when the spec has none. */
  readonly sequence?: ViewerSequence;
}

/** The sequence view's models: the viewer, its source, and the canvas ink it is skinned with (set by the page). */
export interface ViewerSequence {
  readonly viewers: SequenceViewers;
  readonly source: SequenceSource;
  setInk(ink: CanvasInk): void;
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
  /** Raw plugin work as one history step (README §3.2); calls it makes through `api` join the step. */
  run<T>(label: string, work: (plugin: PluginContext, api: ViewerApi) => Promise<T>): Promise<T>;
  readonly managers: ViewerManagers;
  readonly views: PrismViewerViews;
  /** Read-only states over the plugin (state.md §4): what views, view models and the agent's account read. */
  readonly sources: ViewerSources;
  readonly state: {
    readonly ready: ReadonlyState<boolean>;
    readonly busy: ReadonlyState<boolean>;
    /** Nothing loaded. */
    readonly empty: ReadonlyState<boolean>;
  };
  readonly events: {
    readonly hover: Observable<ViewerPointerEvent>;
    readonly click: Observable<ViewerPointerEvent>;
    readonly error: Observable<ViewerError>;
  };
  /** The extensions: which are loaded, each lazy one's status, and loading one (README §3.6). */
  readonly extensions: Pick<ExtensionHost, 'loaded' | 'status' | 'events' | 'load' | 'ensureFor'>;
  /** The session's record (journal.md): what the person and the agent did, in order; in memory until kept. */
  readonly record: SessionRecord;
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
  /** Registers commands with the spec's keybindings applied (README §3.3: keybindings name commands). */
  registerCommands(defs: readonly CommandDefinition<ViewerCommandContext>[]): () => void;
  /** The spec's overrides, for the shell's own parts. */
  readonly overrides: Overrides;
  /** The page's shell takes the extensions' views. */
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

    const targetFormat = createTargetFormat(plugin);
    const history = new HistoryManager<ViewerCheckpoint>(createCheckpointPort(plugin));
    const notifications = new NotificationsManager();
    const { viewState, navigation } = layoutManagers(spec, options);
    const transactions = createUiTransactionCoordinator();
    const managers: ViewerManagers = {
      commands: new CommandManager<ViewerCommandContext>({
        context: commandContext.asReadonly,
        formats: { get: (name) => (name === 'target' ? targetFormat : undefined) },
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

    const targets = createTargetResolver(plugin);
    const modules = new Set(o.modules);
    // Operations that load what an extension owns ask the extensions first; set once they are added.
    const extensionPort: ExtensionPort = { ensureFor: async () => {} };
    const host = createViewerApi({
      plugin,
      sources,
      targets,
      modules: { has: (k) => modules.has(k) },
      managers,
      extensions: extensionPort,
    });
    disposers.push(() => host.dispose());

    const keys = new Map((spec.ui?.keybindings ?? []).map((k) => [k.command, k.key]));
    // Every command goes through here: the core's, the shell's, extensions' and stubs' (overrides.ts).
    const overrides = new Overrides(spec.overrides);
    const registerCommands = (defs: readonly CommandDefinition<ViewerCommandContext>[]) =>
      managers.commands.register(
        overrides.commands(defs).map((d) => {
          if (!keys.has(d.id)) return d;
          const key = keys.get(d.id);
          const { shortcut: _, ...rest } = d;
          return key ? { ...rest, shortcut: key } : rest;
        }),
      );

    // The extensions, in order (viewer-spec.md §3): eager ones activated, lazy ones as stubs.
    const person = host.api.as({ kind: 'person' });
    const viewScope = new SourceScope();
    disposers.push(() => viewScope.dispose());
    const scene = createSceneSource(viewScope, plugin, person);
    const presets = createPresetsModel(plugin, person);
    const componentOptions = createComponentOptionsForm(viewScope, plugin, person);
    const quickStyles = createQuickStylesModel(person);
    const homeActions = createHomeActions(viewScope, plugin, person);
    const sequence =
      spec.sequence === false ? undefined : createSequence(viewScope, plugin, host, o.modules, spec, overrides);
    if (sequence) disposers.push(() => sequence.viewers.dispose());
    disposers.push(registerCommands(coreCommands(person)));
    const extensions = new ExtensionHost({
      host,
      person,
      registerCommands,
      openCommand: (id) => void managers.commands.open(id),
      annotations: host.api.annotations,
      ...(sequence ? { addSequenceBehavior: (b) => sequence.viewers.get().addBehavior(b) } : {}),
      applyPlugin: (parts) => applyExtensionPlugin(plugin, parts),
      overrides,
      task: (label) => {
        const id = `extension:${label}`;
        notifications.putTask({ id, label });
        return () => notifications.removeTask(id);
      },
    });
    disposers.push(() => extensions.dispose());
    // One record per viewer session (journal.md §1), in memory; History reads it, the person may keep it, in the spec's
    // rows when it gives them (`persist.rows`, viewer-spec.md §2.6; typed loosely there, a `FiberRows` here).
    const keepIn = spec.persist?.rows as FiberRows | undefined;
    const record = await createSessionRecord({ api: host.api, undone: history.undone, ...(keepIn ? { keepIn } : {}) });
    disposers.push(() => record.dispose());
    // What is kept and where, with Forget (journal.md §7): the viewer's own sources; the Storage view lists them.
    // The origin that holds it, where there is one (a file: page's is 'null').
    const origin = globalThis.location?.origin;
    const storage = createStorageModel(origin && origin !== 'null' ? { origin } : {});
    const own = viewerStorageSources({ ...(spec.persist ? { persist: spec.persist } : {}), record });
    for (const s of own.sources) disposers.push(storage.register(s));
    disposers.push(() => {
      own.dispose();
      storage.dispose();
    });
    // The data palette the options ask for (README §3.7): before anything is loaded, so no step.
    if (options.canvas.palette !== 'default') await setDataPalette(plugin, options.canvas.palette);
    // The light restore (journal.md §6), written after each step when the spec gives it a store.
    const restore = spec.persist?.restore
      ? createLightRestore(viewScope, plugin, history, spec.persist.restore, {
          notify: (message) => notifications.post(message, 'info'),
        })
      : undefined;
    // Features behind an option (viewer-spec.md §3.1: remote states) are left out unless it is on.
    const gated: Readonly<Record<string, boolean>> = { 'remote-states': options.ui.showRemoteStates };
    await extensions.add((spec.extensions ?? []).filter((e) => gated[e.id] !== false));
    extensionPort.ensureFor = (need) => extensions.ensureFor(need);
    // An override of a command nobody registers is a mistake (viewer-spec.md §4); view parts are checked by spec-check.
    const unknown = overrides
      .unused()
      .filter((k) => !isPartId(k) && !ShellCommandIds.includes(k) && !managers.commands.get(k));
    if (unknown.length) {
      throw new PrismViewerSpecError(unknown.map((k) => `override '${k}': there is no command '${k}'`));
    }

    // What the person sees in the status column: results of the person's and the URL's runs, failures, the plugin's
    // errors, running tasks, undo turning off for a large scene.
    const error = new Subject<ViewerError>();
    const subs = [
      host.runner.finished.subscribe((r) => {
        const result = r.result;
        if (!result || r.nested || r.group) return;
        const def = host.runner.definitions.value.get(r.operation);
        if (!result.ok) {
          if (result.code === 'aborted' || result.code === 'confirm') return;
          if (r.origin.kind === 'person' || r.origin.kind === 'url' || r.origin.kind === 'link') {
            notifications.post(result.hint ? `${result.message} ${result.hint}` : result.message, 'warning');
          }
          error.next({ message: result.message, source: 'api' });
        } else if (
          (r.origin.kind === 'person' || r.origin.kind === 'url') &&
          def &&
          !def.read &&
          // View operations (hover's highlight, focus) say nothing; a lazy stub's run is the real one's result.
          (def.history !== false || isStubDefinition(def))
        ) {
          notifications.post(result.did, 'info');
        } else if (r.origin.kind === 'person' && def && !def.read && /^(select|focus|measure)\./.test(r.operation)) {
          // A selection, a focus or a measurement by target is how a person without sight picks (README §3.7): its
          // result is said (the live region), with no message on show.
          notifications.announce(result.did);
        }
      }),
      plugin.events.log.subscribe((e) => {
        if (e.type !== 'error') return;
        notifications.post(String(e.message), 'error');
        error.next({ message: String(e.message), source: 'plugin' });
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
    let disposed = false;
    const core: ViewerCore = {
      plugin,
      api: host.api,
      run: (label, work) =>
        host.runner.opaque(label, { kind: 'script' }, (context) =>
          work(plugin, host.bind({ origin: { kind: 'script' }, within: context })),
        ),
      managers,
      views: { scene, presets, componentOptions, quickStyles, homeActions, ...(sequence ? { sequence } : {}) },
      sources,
      state: { ready: ready.asReadonly, busy: sources.busy, empty },
      events: { hover: sources.pointer.hover, click: sources.pointer.click, error: error.asObservable() },
      extensions,
      record,
      ...(restore ? { restore } : {}),
      storage,
      spec,
      options,
      host,
      ready$: ready,
      registerCommands,
      overrides,
      attachExtensionUI: (port) => extensions.attachUI(port),
      mountCanvas: async (el) => {
        const ok = await mountCanvas(plugin, el);
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

/** Frames for the sequence view: the page's, or timers headless. */
function sequenceFrames(page: boolean): Scheduler {
  return rafScheduler({
    requestAnimationFrame: page
      ? (cb) => requestAnimationFrame(cb)
      : (cb) => setTimeout(() => cb(performance.now()), 16) as unknown as number,
    cancelAnimationFrame: page ? (id) => cancelAnimationFrame(id) : (id) => clearTimeout(id),
    setTimeout: (cb, ms) => setTimeout(cb, ms) as unknown as number,
    clearTimeout: (id) => clearTimeout(id),
    now: () => performance.now(),
  });
}

/**
 * The sequence (sequence.md, plan 2.4): its viewer over the structure the source follows, its area in the API, and its
 * state in snapshots (`prism-viewer.sequence`, sequence.md §7.1).
 */
function createSequence(
  scope: SourceScope,
  plugin: PluginContext,
  host: ViewerApiHost,
  modules: readonly ApiModuleKey[],
  spec: PrismViewerSpec,
  overrides: Overrides,
): ViewerSequence {
  const link = host.api.as({ kind: 'link' });
  const source = createSequenceSource(scope, plugin, link);
  const ink = observableState<CanvasInk>({ theme: 'light', color: '#ffffff' });
  const viewers = createSequenceViewers({
    ink: ink.asReadonly,
    structure: source.structure,
    actions: source.actions,
    scheduler: sequenceFrames(modules.includes('page')),
    // The viewer's own behaviors, overrides applied; the extensions' attach as they activate (sequence.md §3).
    behaviors: overrides.parts(spec.sequence ? spec.sequence.behaviors : []),
  });
  const viewer = viewers.get();
  const initial = spec.sequence ? spec.sequence.style : undefined;
  if (initial === 'text') viewer.setStyle('compact');
  scope.add(host.addArea('sequence', (ctx) => new SequenceApi(ctx, viewer, source)));
  scope.add(
    plugin.state.registerCustomState('prism-viewer.sequence', {
      get: () => ({
        version: 1,
        style: viewer.style.value,
        colors: viewer.theme.value.letters,
        rail: viewer.rail.value,
        structure: source.followed.value,
      }),
      set: (value) => {
        const v = value as
          | { version: number; style: SequenceStyle; colors: string | null; rail: TracksRail; structure?: string }
          | undefined;
        if (!v || v.version !== 1) {
          viewer.setStyle(initial === 'text' ? 'compact' : 'tracks');
          viewer.setTheme(DefaultSequenceTheme);
          source.setFollowed(undefined);
          return;
        }
        viewer.setStyle(v.style);
        viewer.setTheme({ letters: v.colors });
        if (v.rail) viewer.setRail(v.rail);
        // Undo's checkpoints leave out what the sequence follows (view state); a session has it.
        if ('structure' in v) source.setFollowed(v.structure);
      },
    }),
  );
  return { viewers, source, setInk: (i) => ink.set(i) };
}
