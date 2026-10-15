/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { PluginContext } from '@molstar/plugin/context';
import type { ViewerSources } from '../sources/viewer-sources.js';
import {
  type ApiArea,
  type ApiAreaContext,
  type ApiModules,
  ApiResult,
  BuiltInAreaFeatures,
  bindOperation,
  type Changed,
  type ExtensionPort,
  type GroupOptions,
  type GroupResult,
  isOperationFn,
  type JournalPort,
  NoJournal,
  type OperationFn,
  type RunInput,
  type RunOptions,
  type RunOrigin,
  type ViewerApi,
} from './area.js';
import { CameraApi } from './camera.js';
import { CanvasApi } from './canvas.js';
import { FocusApi } from './focus.js';
import { HighlightApi } from './highlight.js';
import { HistoryApi } from './history.js';
import { LoadApi } from './load.js';
import type { ViewerManagers } from './managers.js';
import { ApiRunner, messageOf, type RunRecord } from './run.js';
import { SelectApi } from './select.js';
import { SessionApi } from './session.js';
import { SettingsApi } from './settings.js';
import type { TargetResolver } from './targets.js';
import { ViewApi } from './view.js';

// The viewer's API object (api.md §2–§6). Each area is made once; the API a caller gets is a facade bound to how its
// calls run: `viewer.api` as `script`, `as(origin)` for each surface (view models and commands `person`, tools `agent`
// with the turn, the loaders feature `url`, the sequence `link`), `group` for its calls, `ctx.api` for an operation's
// own calls (part of its step). Facades are cheap: an area's operations are rebound on first use.

/** Makes an area, once per viewer. */
export type AreaFactory = (ctx: ApiAreaContext) => ApiArea;

/** The core areas (api.md §2), as they are built; the built-in features bring theirs (`BuiltInApis`). */
export const CoreAreas: Readonly<Record<string, AreaFactory>> = {
  canvas: (c) => new CanvasApi(c),
  camera: (c) => new CameraApi(c),
  load: (c) => new LoadApi(c),
  session: (c) => new SessionApi(c),
  select: (c) => new SelectApi(c),
  focus: (c) => new FocusApi(c),
  highlight: (c) => new HighlightApi(c),
  history: (c) => new HistoryApi(c),
  view: (c) => new ViewApi(c),
  settings: (c) => new SettingsApi(c),
};

export interface ViewerApiDeps {
  readonly plugin: PluginContext;
  readonly sources: ViewerSources;
  readonly targets: TargetResolver;
  readonly modules: ApiModules;
  readonly managers: ViewerManagers;
  readonly extensions: ExtensionPort;
  /** The record's part of session files; none: `NoJournal`. */
  readonly journal?: JournalPort;
}

export interface ViewerApiHost {
  /** The API as a script calls it (origin `script`). */
  readonly api: ViewerApi;
  readonly runner: ApiRunner;
  /** The context areas are made with. */
  readonly areaContext: ApiAreaContext;
  /** The API bound to how calls run (an origin, a run or group in progress). */
  bind(input: RunInput): ViewerApi;
  /** Adds an area (an extension's); returns its removal. Its id must be new. */
  addArea(id: string, factory: AreaFactory): () => void;
  /** Replaces an area in place (a lazy extension's stub by the real one); returns its removal. */
  replaceArea(id: string, factory: AreaFactory): () => void;
  dispose(): void;
}

const ApiMembers = new Set(['run', 'as', 'group', 'area', 'describe', 'definitions', 'running', 'events']);

export function createViewerApi(deps: ViewerApiDeps): ViewerApiHost {
  const runner = new ApiRunner();
  const areaContext: ApiAreaContext = {
    plugin: deps.plugin,
    sources: deps.sources,
    targets: deps.targets,
    modules: deps.modules,
    managers: deps.managers,
    extensions: deps.extensions,
    journal: deps.journal ?? NoJournal,
    operations: runner,
  };
  const areas = new Map<string, ApiArea>();
  const cache = new Map<string, ViewerApi>();

  const make = (id: string, factory: AreaFactory) => {
    const area = factory(areaContext);
    if (area.id !== id) throw new Error(`The area made for '${id}' says it is '${area.id}'`);
    return area;
  };
  for (const [id, factory] of Object.entries(CoreAreas)) areas.set(id, make(id, factory));

  /** An area seen through a binding: its operations rebound, everything else its own. */
  const facadeOf = (area: ApiArea, input: RunInput): ApiArea => {
    const bound = new Map<string, OperationFn<any, any>>();
    return new Proxy(area, {
      get(target, prop, receiver) {
        const value = Reflect.get(target, prop, target);
        if (typeof prop === 'string' && isOperationFn(value)) {
          let fn = bound.get(prop);
          if (!fn || fn.definition !== value.definition) {
            fn = bindOperation(runner, value.definition, input);
            bound.set(prop, fn);
          }
          return fn;
        }
        return typeof value === 'function' && prop !== 'constructor' ? value.bind(target) : value;
      },
    });
  };

  const bind = (input: RunInput): ViewerApi => {
    // Plain origins are kept (one per kind); turns, groups and runs in progress are made each time (they are cheap).
    const key = !input.within && !input.group && input.origin.kind !== 'agent' ? input.origin.kind : undefined;
    if (key) {
      const cached = cache.get(key);
      if (cached) return cached;
    }
    const facades = new Map<ApiArea, ApiArea>();
    const facade = (area: ApiArea) => {
      let f = facades.get(area);
      if (!f) {
        f = facadeOf(area, input);
        facades.set(area, f);
      }
      return f;
    };
    const members = {
      run: (id: string, values?: Readonly<Record<string, unknown>>, options?: RunOptions) => {
        const def = runner.definitions.value.get(id);
        if (!def) {
          return Promise.resolve(
            ApiResult.fail(
              'not-found',
              `There is no operation '${id}'.`,
              'List the operations with the tools or the command menu.',
            ),
          );
        }
        return runner.run(def, values ?? {}, { ...input, ...options });
      },
      as: (origin: RunOrigin) => bind({ ...input, origin }),
      group: <T>(label: string, work: (api: ViewerApi) => Promise<T>, options: GroupOptions = {}) =>
        runGroup(input, label, work, options),
      area: (id: string) => {
        const area = areas.get(id);
        return area ? facade(area) : undefined;
      },
      describe: (id: string) => runner.describe(id),
      definitions: runner.definitions,
      running: runner.running,
      events: { started: runner.started, finished: runner.finished },
    };
    const api = new Proxy(members, {
      get(target, prop) {
        if (typeof prop === 'string' && !ApiMembers.has(prop)) {
          const area = areas.get(prop);
          if (area) return facade(area);
          // A built-in feature's area this viewer's spec left out: say which feature brings it (`BuiltInApis`).
          const feature = BuiltInAreaFeatures[prop as keyof typeof BuiltInAreaFeatures];
          if (feature) {
            throw new Error(
              `The '${prop}' area is not in this viewer: add ${feature} to its spec's extensions (DefaultPrismViewerSpec has it).`,
            );
          }
        }
        return Reflect.get(target, prop);
      },
      has(target, prop) {
        return (typeof prop === 'string' && areas.has(prop)) || Reflect.has(target, prop);
      },
      ownKeys(target) {
        return [...Reflect.ownKeys(target), ...areas.keys()];
      },
      getOwnPropertyDescriptor(target, prop) {
        if (typeof prop === 'string' && areas.has(prop)) {
          return { configurable: true, enumerable: true, value: facade(areas.get(prop)!) };
        }
        return Reflect.getOwnPropertyDescriptor(target, prop);
      },
    }) as unknown as ViewerApi;
    if (key) cache.set(key, api);
    return api;
  };

  const runGroup = async <T>(
    input: RunInput,
    label: string,
    work: (api: ViewerApi) => Promise<T>,
    options: GroupOptions,
  ): Promise<GroupResult<T>> => {
    const group = runner.newGroup(label);
    const history = deps.managers.history;
    const within = history.isActive(input.within) ? input.within : undefined;
    const failure = (f: Extract<ApiResult, { ok: false }>, restored: boolean): GroupResult<T> => ({
      ...f,
      hint: f.hint ?? (restored ? 'Nothing was changed: the whole group was undone.' : undefined),
      runs: group.runs,
    });
    try {
      const outcome = await history.group(
        label,
        input.origin,
        (context) =>
          work(bind({ ...input, within: context, group, ...(options.signal ? { signal: options.signal } : {}) })),
        {
          failed: (v) => ApiResult.isFailure(v),
          onFailure: options.onFailure ?? 'rollback',
          ...(within ? { within } : {}),
        },
      );
      if (outcome.failed) {
        const first = ApiResult.isFailure(outcome.value)
          ? outcome.value
          : ((group.runs.find((r) => r.result && !r.result.ok)?.result as Extract<ApiResult, { ok: false }>) ??
            ApiResult.fail('failed', `${label} failed.`));
        return failure(first, outcome.restored);
      }
      return {
        ok: true,
        did: label,
        changed: unionOf(group.runs),
        value: outcome.value,
        runs: group.runs,
        ...(outcome.runs.length ? {} : { unchanged: true as const }),
      };
    } catch (e) {
      return failure(ApiResult.fail('failed', `${label} failed: ${messageOf(e)}`), true);
    }
  };

  runner.connect({ ...deps, bind }, areaContext);
  const scriptApi = bind({ origin: { kind: 'script' } });

  const remove = (id: string, area: ApiArea) => {
    if (areas.get(id) !== area) return;
    for (const op of area.operations()) runner.undefine(op.definition.id);
    area.dispose?.();
    areas.delete(id);
  };

  return {
    api: scriptApi,
    runner,
    areaContext,
    bind,
    addArea(id, factory) {
      if (areas.has(id)) throw new Error(`An API area '${id}' exists already`);
      const area = make(id, factory);
      areas.set(id, area);
      return () => remove(id, area);
    },
    replaceArea(id, factory) {
      const old = areas.get(id);
      if (old) remove(id, old);
      const area = make(id, factory);
      areas.set(id, area);
      return () => remove(id, area);
    },
    dispose() {
      runner.dispose();
      for (const [id, area] of [...areas]) remove(id, area);
      cache.clear();
    },
  };
}

/** What a group's calls changed, together. */
function unionOf(runs: readonly RunRecord[]): Changed {
  const nodes = new Set<string>();
  let selection = false;
  let view = false;
  let sequence = false;
  for (const r of runs) {
    if (!r.result?.ok) continue;
    const c = r.result.changed;
    for (const n of c.nodes ?? []) nodes.add(n);
    selection ||= !!c.selection;
    view ||= !!c.view;
    sequence ||= !!c.sequence;
  }
  return {
    ...(nodes.size ? { nodes: [...nodes] } : {}),
    ...(selection ? { selection } : {}),
    ...(view ? { view } : {}),
    ...(sequence ? { sequence } : {}),
  };
}

export * from './area.js';
export type { RunRecord } from './run.js';
export { ApiRunner, checkValues } from './run.js';
export type { ViewerCommandContext, ViewerManagers } from './managers.js';
export type { ResolvedPart, ResolvedTarget, ResolveOptions, TargetResolver } from './targets.js';
