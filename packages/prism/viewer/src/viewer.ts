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
  createUiTransactionCoordinator,
  HistoryManager,
  NavigationManager,
  NotificationsManager,
  observableState,
  type ReadonlyState,
  ViewStateManager,
} from '@molstar/prism-ui';
import { type Observable, Subject } from 'rxjs';
import type { ApiModuleKey, ViewerApi } from './api/area.js';
import { createViewerApi, type ViewerApiHost } from './api/index.js';
import type { ViewerCommandContext, ViewerManagers } from './api/managers.js';
import { NoTargets, type TargetResolver } from './api/targets.js';
import { coreCommands } from './commands/core.js';
import { type PrismViewerOptionValues, resolveOptions } from './options.js';
import { createCheckpointPort, type ViewerCheckpoint, watchOutsideChanges } from './sources/checkpoint.js';
import { composePluginSpec } from './sources/compose-plugin.js';
import { createViewerSources, type ViewerSources } from './sources/viewer-sources.js';
import { isLazyExtension, type PrismViewerSpec, type ViewerExtension } from './spec.js';
import { checkSpec } from './spec-check.js';

// The viewer (viewer-spec.md §6): created from a `PrismViewerSpec`, in a fixed order: validate the spec, compose the
// plugin spec, create and `init()` the plugin, create the sources, managers and API, activate the extensions (lazy ones
// as stubs); a page then mounts the canvas and renders the shell (`createPrismViewer`), headless has neither
// (`createPrismViewerHeadless`). `ready` is set last. `dispose()` undoes it in reverse.

/** The root view models (README §3.3): the shell's on a page; the scene, the selection and the rest as they are built. */
export interface PrismViewerViews {
  /** The shell (docks, tabs, cards); undefined headless. */
  readonly shell?: import('@molstar/prism-ui').ShellModel;
}

/** What is under the pointer, for embedding hosts (v5's `structureInteractivity`); targets with plan phase 2.3. */
export interface ViewerPointerEvent {
  /** The plugin's label of it, as lines; empty: nothing. */
  readonly lines: readonly string[];
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
  /** Raw plugin work as one history step (README §3.2). */
  run<T>(label: string, work: (plugin: PluginContext) => Promise<T>): Promise<T>;
  readonly managers: ViewerManagers;
  readonly views: PrismViewerViews;
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
  readonly spec: PrismViewerSpec;
  readonly options: PrismViewerOptionValues;
  dispose(): Promise<void>;
}

export interface ViewerCoreOptions {
  readonly spec: PrismViewerSpec;
  /** Makes the plugin from the composed spec (a `HeadlessPluginContext` when images are wanted headless). */
  readonly createPlugin: (spec: PluginSpec) => PluginContext;
  /** What the host offers (api.md §4). */
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
}

export async function createViewerCore(o: ViewerCoreOptions): Promise<ViewerCore> {
  const { spec } = o;
  checkSpec(spec);
  const options = resolveOptions(spec.options);
  const eager = (spec.extensions ?? []).filter((e): e is ViewerExtension => !isLazyExtension(e));
  const pluginSpec = composePluginSpec({ base: spec.plugin, extensions: eager, options, finalize: spec.finalize });

  const plugin = o.createPlugin(pluginSpec);
  await plugin.init();

  const disposers: (() => void | Promise<void>)[] = [];
  const sources = createViewerSources(plugin);
  disposers.push(() => sources.dispose());

  const empty = sources.empty;
  const commandContext = observableState<ViewerCommandContext>({ empty: empty.value });
  const emptySub = empty.subscribe((e) => commandContext.set({ empty: e }));
  disposers.push(() => emptySub.unsubscribe());

  const history = new HistoryManager<ViewerCheckpoint>(createCheckpointPort(plugin));
  const notifications = new NotificationsManager();
  const viewState = new ViewStateManager();
  const navigation = new NavigationManager();
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

  const targets: TargetResolver = NoTargets;
  const modules = new Set(o.modules);
  const host = createViewerApi({ plugin, sources, targets, modules: { has: (k) => modules.has(k) }, managers });
  disposers.push(() => host.dispose());

  // The extensions, in order: their API areas and commands (their views with the shell, lazy ones as stubs: plan 3.1).
  const person = host.api.as({ kind: 'person' });
  for (const e of eager) {
    if (e.api) disposers.push(host.addArea(e.id, e.api));
  }
  disposers.push(managers.commands.register(coreCommands(person)));
  for (const e of eager) {
    const commands = e.commands?.(person);
    if (commands?.length) disposers.push(managers.commands.register(commands));
  }

  // What the person sees in the status column: results of the person's and the URL's runs, failures, the plugin's
  // errors, running tasks, undo turning off for a large scene.
  const error = new Subject<ViewerError>();
  const subs = [
    host.runner.finished.subscribe((r) => {
      const result = r.result;
      if (!result) return;
      const def = host.runner.definitions.get(r.operation);
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

  const hover = new Subject<ViewerPointerEvent>();
  const click = new Subject<ViewerPointerEvent>();
  disposers.push(() => {
    hover.complete();
    click.complete();
  });

  const ready = observableState(false);
  let disposed = false;
  const core: ViewerCore = {
    plugin,
    api: host.api,
    run: (label, work) => host.runner.opaque(label, { kind: 'script' }, () => work(plugin)),
    managers,
    views: {},
    sources,
    state: { ready: ready.asReadonly, busy: sources.busy, empty },
    events: { hover: hover.asObservable(), click: click.asObservable(), error: error.asObservable() },
    spec,
    options,
    host,
    ready$: ready,
    extensions: eager,
    onDispose: (fn) => disposers.push(fn),
    async dispose() {
      if (disposed) return;
      disposed = true;
      ready.set(false);
      for (const d of disposers.reverse()) {
        try {
          await d();
        } catch (e) {
          console.error(e);
        }
      }
      plugin.dispose();
    },
  };
  return core;
}
