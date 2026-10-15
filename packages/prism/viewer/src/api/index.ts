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
  type RunOptions,
  type RunOrigin,
  type ViewerApi,
} from './area.js';
import { CameraApi } from './camera.js';
import { CanvasApi } from './canvas.js';
import { HistoryApi } from './history.js';
import { LoadApi } from './load.js';
import type { ViewerManagers } from './managers.js';
import { ApiRunner } from './run.js';
import { SessionApi } from './session.js';
import { StructureApi } from './structure.js';
import type { TargetResolver } from './targets.js';
import { ViewApi } from './view.js';

// The viewer's API object (api.md §2–§6): the core areas, extensions' areas added by the composition, bound by origin.
// `viewer.api` runs as `script`; each surface binds its own (`as`): view models and commands `person`, tools `agent`
// with the turn, the loaders feature `url`, the sequence `link`.

/** Makes an area bound to an origin (and an `api.group`'s token): one per origin, made when first bound. */
export type AreaFactory = (ctx: ApiAreaContext, origin: RunOrigin, group?: object) => ApiArea;

/** The core areas (api.md §2), as they are built. */
export const CoreAreas: Readonly<Record<string, AreaFactory>> = {
  canvas: (c, o, g) => new CanvasApi(c, o, g),
  camera: (c, o, g) => new CameraApi(c, o, g),
  load: (c, o, g) => new LoadApi(c, o, g),
  session: (c, o, g) => new SessionApi(c, o, g),
  structure: (c, o, g) => new StructureApi(c, o, g),
  history: (c, o, g) => new HistoryApi(c, o, g),
  view: (c, o, g) => new ViewApi(c, o, g),
};

export interface ViewerApiDeps {
  readonly plugin: PluginContext;
  readonly sources: ViewerSources;
  readonly targets: TargetResolver;
  readonly modules: ApiModules;
  readonly managers: ViewerManagers;
}

export interface ViewerApiHost {
  /** The API as a script calls it (origin `script`). */
  readonly api: ViewerApi;
  readonly runner: ApiRunner;
  /** Adds an area (an extension's); returns its removal. */
  addArea(id: string, factory: AreaFactory): () => void;
  dispose(): void;
}

function originKey(origin: RunOrigin): string {
  return origin.kind === 'agent' ? `agent:${origin.turn}` : origin.kind;
}

export function createViewerApi(deps: ViewerApiDeps): ViewerApiHost {
  const runner = new ApiRunner();
  const areaCtx: ApiAreaContext = { ...deps, runner };
  const factories = new Map<string, AreaFactory>(Object.entries(CoreAreas));
  const bound = new Map<string, { origin: RunOrigin; api: ViewerApi }>();
  const make = (factory: AreaFactory, origin: RunOrigin, group?: object) => factory(areaCtx, origin, group);

  const bind = (origin: RunOrigin, group?: object): ViewerApi => {
    const key = originKey(origin);
    if (!group) {
      const cached = bound.get(key);
      if (cached) return cached.api;
    }
    const areas = new Map<string, ApiArea>();
    for (const [id, factory] of factories) areas.set(id, make(factory, origin, group));
    const api = {
      run: (id: string, values?: Readonly<Record<string, unknown>>, options?: RunOptions) => {
        const def = runner.definitions.get(id);
        if (!def) {
          return Promise.resolve(
            ApiResult.fail(
              'not-found',
              `There is no operation '${id}'.`,
              'List the operations with the tools or the command menu.',
            ),
          );
        }
        return runner.run(def, values ?? {}, { origin, ...(group ? { group } : {}), ...options });
      },
      as: (other: RunOrigin) => bind(other),
      group: async <T>(label: string, work: (api: ViewerApi) => Promise<T>): Promise<ApiResult<T>> => {
        const token = {};
        runner.nameGroup(token, label);
        try {
          const outcome = await deps.managers.history.group(label, origin, () => work(bind(origin, token)), {
            failed: (v) => isFailure(v),
          });
          if (isFailure(outcome.value)) {
            const f = outcome.value;
            return { ...f, hint: f.hint ?? 'Nothing was changed: the whole group was undone.' };
          }
          return ApiResult.ok(label, {}, { value: outcome.value });
        } catch (e) {
          return ApiResult.fail(
            'failed',
            `${label} failed: ${e instanceof Error ? e.message : String(e)}`,
            'Nothing was changed: the whole group was undone.',
          );
        }
      },
      area: (id: string) => areas.get(id),
      running: runner.running,
      events: { started: runner.started, finished: runner.finished },
    } as Record<string, unknown>;
    for (const [id, area] of areas) {
      Object.defineProperty(api, id, { value: area, enumerable: true, configurable: true });
    }
    const result = api as unknown as ViewerApi;
    if (!group) bound.set(key, { origin, api: result });
    return result;
  };

  runner.connect({ ...deps, bind });
  const scriptApi = bind({ kind: 'script' });

  return {
    api: scriptApi,
    runner,
    addArea(id, factory) {
      if (factories.has(id)) throw new Error(`An API area '${id}' exists already`);
      factories.set(id, factory);
      for (const { origin, api } of bound.values()) {
        Object.defineProperty(api, id, { value: make(factory, origin), enumerable: true, configurable: true });
      }
      return () => {
        if (factories.get(id) !== factory) return;
        factories.delete(id);
        for (const { api } of bound.values()) {
          const area = (api as unknown as Record<string, ApiArea | undefined>)[id];
          for (const op of area?.operations() ?? []) runner.undefine(op.definition.id);
          delete (api as unknown as Record<string, unknown>)[id];
        }
      };
    },
    dispose() {
      runner.dispose();
      bound.clear();
    },
  };
}

function isFailure(v: unknown): v is Extract<ApiResult, { ok: false }> {
  return !!v && typeof v === 'object' && (v as { ok?: unknown }).ok === false && 'code' in (v as object);
}

export * from './area.js';
export type { RunRecord } from './run.js';
export { ApiRunner, checkValues } from './run.js';
export type { ViewerCommandContext, ViewerManagers } from './managers.js';
export type { ResolvedPart, ResolvedTarget, ResolveOptions, TargetResolver } from './targets.js';
