/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { PluginContext } from '@molstar/plugin/context';
import type { PluginSpec } from '@molstar/plugin/spec';
import {
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
import type { ApiModuleKey, ViewerApi } from './api/area.js';
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
import { SourceScope } from './sources/state.js';
import { createComponentOptionsForm, createPresetsModel, createQuickStylesModel } from './sources/structure-tab.js';
import { createTargetResolver } from './sources/targets.js';
import { createViewerSources, type ViewerSources } from './sources/viewer-sources.js';
import { isLazyExtension, type PrismViewerSpec, type ViewerExtension } from './spec.js';
import { checkSpec } from './spec-check.js';
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
  /** The eager extensions, in order. */
  readonly extensions: readonly ViewerExtension[];
  /** Registers commands with the spec's keybindings applied (README §3.3: keybindings name commands). */
  registerCommands(defs: readonly CommandDefinition<ViewerCommandContext>[]): () => void;
  /** Mounts the canvas in a page's element; false (and the `canvas` module dropped) when it cannot start. */
  mountCanvas(host: HTMLElement): Promise<boolean>;
}

/** The toolbar items `options.ui.buttons` hides, by option (viewer-spec.md §2.2). */
const ButtonIds: Readonly<Record<string, string>> = {
  reset: 'reset',
  screenshot: 'screenshot',
  expand: 'expand',
  fullscreen: 'fullscreen',
  settings: 'settings',
  selectionMode: 'selection-mode',
  animation: 'animation',
  trajectory: 'trajectory',
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

    const history = new HistoryManager<ViewerCheckpoint>(createCheckpointPort(plugin));
    const notifications = new NotificationsManager();
    const { viewState, navigation } = layoutManagers(spec, options);
    const transactions = createUiTransactionCoordinator();
    const managers: ViewerManagers = {
      commands: new CommandManager<ViewerCommandContext>({
        context: commandContext.asReadonly,
        formats: { get: () => undefined },
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
    const host = createViewerApi({ plugin, sources, targets, modules: { has: (k) => modules.has(k) }, managers });
    disposers.push(() => host.dispose());

    const keys = new Map((spec.ui?.keybindings ?? []).map((k) => [k.command, k.key]));
    const registerCommands = (defs: readonly CommandDefinition<ViewerCommandContext>[]) =>
      managers.commands.register(
        defs.map((d) => {
          if (!keys.has(d.id)) return d;
          const key = keys.get(d.id);
          const { shortcut: _, ...rest } = d;
          return key ? { ...rest, shortcut: key } : rest;
        }),
      );

    // The extensions, in order: their API areas, contributed operations and commands (their views with the shell,
    // lazy ones as stubs: plan 3.1).
    const person = host.api.as({ kind: 'person' });
    const viewScope = new SourceScope();
    disposers.push(() => viewScope.dispose());
    const scene = createSceneSource(viewScope, plugin, person);
    const presets = createPresetsModel(plugin, person);
    const componentOptions = createComponentOptionsForm(viewScope, plugin, person);
    const quickStyles = createQuickStylesModel(person);
    const homeActions = createHomeActions(viewScope, plugin, person);
    for (const e of eager) {
      for (const [id, factory] of Object.entries(e.api ?? {})) disposers.push(host.addArea(id, factory));
      for (const c of e.contributes?.(host.areaContext) ?? []) {
        if (!e.builtIn)
          throw new Error(
            `Extension '${e.id}' adds '${c.operation.id}' to a core area: only the viewer's features may`,
          );
        const area = host.api.area(c.area);
        if (!area) throw new Error(`Extension '${e.id}' adds to an area that does not exist: '${c.area}'`);
        disposers.push(area.contribute(c.name, c.operation));
      }
    }
    for (const e of eager) {
      for (const a of e.annotations ?? []) disposers.push(host.api.annotations.register(a));
    }
    disposers.push(registerCommands(coreCommands(person)));
    for (const e of eager) {
      const commands = e.commands?.(person);
      if (commands?.length) disposers.push(registerCommands(commands));
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
        } else if ((r.origin.kind === 'person' || r.origin.kind === 'url') && def && !def.read) {
          notifications.post(result.did, 'info');
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
      views: { scene, presets, componentOptions, quickStyles, homeActions },
      sources,
      state: { ready: ready.asReadonly, busy: sources.busy, empty },
      events: { hover: sources.pointer.hover, click: sources.pointer.click, error: error.asObservable() },
      spec,
      options,
      host,
      ready$: ready,
      extensions: eager,
      registerCommands,
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
    docks[id] = {
      ...(tabs ? { tabs } : {}),
      ...(size !== undefined ? { size } : {}),
      ...(floating !== undefined ? { floating } : {}),
      ...(pinned !== undefined ? { pinned } : {}),
    };
    nav[id] = { ...(active !== undefined ? { active } : {}), ...(collapsed !== undefined ? { collapsed } : {}) };
  }
  // A tab placed by the spec: listed in its dock, by order, so it lands there when it registers.
  const placed = Object.entries(layout.tabs ?? {}).sort((a, b) => (a[1].order ?? 0) - (b[1].order ?? 0));
  for (const [tab, place] of placed) {
    const dock = place.dock ?? 'right';
    const current = docks[dock]?.tabs ?? [];
    docks[dock] = { ...docks[dock], tabs: [...current.filter((t) => t !== tab), tab] };
  }
  const buttons: Record<string, boolean> = {};
  for (const [option, on] of Object.entries(options.ui.buttons)) if (!on) buttons[ButtonIds[option] ?? option] = false;
  const viewState = new ViewStateManager({
    docks,
    ...(layout.mode ? { mode: layout.mode } : {}),
    compactRail: !!layout.compactRail,
    buttons,
  });
  const navigation = new NavigationManager({ docks: nav });
  return { viewState, navigation };
}
