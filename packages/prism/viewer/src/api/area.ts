/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import type { PluginContext } from '@molstar/plugin/context';
import type { ReadonlyState } from '@molstar/prism-ui';
import type { Observable } from 'rxjs';
import type { ApiRunner, RunRecord } from './run.js';
import type { ViewerManagers } from './managers.js';
import type { TargetResolver } from './targets.js';
import type { ViewerSources } from '../sources/viewer-sources.js';

// The viewer's API (api.md, operations.md §1): one object per area, its operations as fields, every change through
// one run path (`ApiRunner`), results for every caller, never a throw.

/** Who runs it (operations.md §1): the history and the record read it. */
export type RunOrigin =
  | { readonly kind: 'person' }
  | { readonly kind: 'agent'; readonly turn: string }
  | { readonly kind: 'script' }
  | { readonly kind: 'url' }
  | { readonly kind: 'link' };

export type RunOriginKind = RunOrigin['kind'];

/** What a reader re-reads after a run. */
export interface Changed {
  readonly nodes?: readonly string[];
  readonly selection?: boolean;
  readonly view?: boolean;
  readonly sequence?: boolean;
}

export type ApiResult<R = unknown> =
  | {
      readonly ok: true;
      /** One line for the person or the model: what happened ("Background set to white."). */
      readonly did: string;
      readonly changed: Changed;
      readonly value?: R;
      readonly warnings?: readonly string[];
      /** What to do next, when there is an obvious next step. */
      readonly next?: string;
    }
  | {
      readonly ok: false;
      /** `invalid`, `blocked`, `unavailable`, `not-found`, `aborted`, `failed`, … */
      readonly code: string;
      readonly message: string;
      /** What to do instead. */
      readonly hint?: string;
    };

export type ApiFailure = Extract<ApiResult, { ok: false }>;

export const ApiResult = {
  ok<R>(
    did: string,
    changed: Changed = {},
    extra: { value?: R; warnings?: string[]; next?: string } = {},
  ): ApiResult<R> {
    return { ok: true, did, changed, ...extra };
  },
  fail(code: string, message: string, hint?: string): ApiFailure {
    return hint ? { ok: false, code, message, hint } : { ok: false, code, message };
  },
};

/** What the host offers; an operation's `modules` are checked against it (api.md §4). */
export type ApiModuleKey = 'canvas' | 'files' | 'network' | 'storage' | 'page';

export interface ApiModules {
  has(key: ApiModuleKey): boolean;
}

/** An optimistic edit (`managers.transactions`): applied to a view model's base until the source reports it. */
export interface Projection<T = unknown> {
  /** Which view model's base it projects onto (its id in the projection registry). */
  readonly model: string;
  apply(base: T): T;
  isAbsorbed(base: T): boolean;
}

/** What an operation's run gets (api.md §4). */
export interface ApiContext {
  /** Only `api/` and `sources/` use it (contexts.md §7). */
  readonly plugin: PluginContext;
  readonly sources: ViewerSources;
  readonly targets: TargetResolver;
  /** Other areas, bound to the same origin (and group). */
  readonly api: ViewerApi;
  readonly origin: RunOrigin;
  readonly modules: ApiModules;
  readonly managers: ViewerManagers;
  readonly signal: AbortSignal;
  /** Asked for another of what exists (decision 30). */
  readonly copy: boolean;
}

/** What creating an area gets: the shared parts, no origin. */
export interface ApiAreaContext {
  readonly plugin: PluginContext;
  readonly sources: ViewerSources;
  readonly targets: TargetResolver;
  readonly modules: ApiModules;
  readonly managers: ViewerManagers;
  readonly runner: ApiRunner;
}

export interface OperationDefinition<P extends PD.Params = PD.Params, R = unknown> {
  /** `<area>.<verb>`: `component.add`, `canvas.background`. */
  readonly id: string;
  /** For people: the command's title, and the history entry when the result says nothing better. */
  readonly title: string;
  /** One definition: the command's form and the tool's schema derive from it. */
  readonly params: P | ((ctx: ApiContext) => P);
  /** Changes nothing: no history step, may run alongside others. */
  readonly read?: true;
  /** Changes view state, not content: no history step (select, focus, camera, highlight, view). */
  readonly history?: false;
  /** A step the history cannot restore (an upload); default true. */
  readonly undoable?: false;
  /** Asks before it runs, from the agent (agent.md §4). */
  readonly approval?: 'network' | 'upload' | 'destructive';
  /** What it needs of the host; where one is missing it says why instead of running. */
  readonly modules?: readonly ApiModuleKey[];
  /** Whether it can run now; the string is the reason. */
  when?(ctx: ApiContext): true | string;
  /**
   * Identity for operations that create things (decision 30): the thing it would create already exists, said as a
   * result line ("1TQN is already loaded."); the run then does nothing, unless `copy` asks for another.
   */
  exists?(values: PD.ValuesFor<P>, ctx: ApiContext): string | undefined;
  /** An optimistic edit, when the result is known before the run ends. */
  project?(values: PD.ValuesFor<P>): Projection;
  /** Who may not call it: `ui-only` / `agent-only` say why there is no tool or no command (operations.md §4). */
  readonly surface?: { readonly kind: 'ui-only' | 'agent-only'; readonly reason: string };
  run(values: PD.ValuesFor<P>, ctx: ApiContext): Promise<ApiResult<R>>;
}

export interface RunOptions {
  readonly signal?: AbortSignal;
  /** Another of what exists (decision 30). */
  readonly copy?: boolean;
}

export interface OperationFn<P extends PD.Params = PD.Params, R = unknown> {
  (values?: Partial<PD.ValuesFor<P>>, options?: RunOptions): Promise<ApiResult<R>>;
  readonly definition: OperationDefinition<P, R>;
}

/** The base of every area: its id, and `operation` to declare an operation as a field (api.md §3). */
export abstract class ApiArea {
  abstract readonly id: string;

  constructor(
    protected readonly ctx: ApiAreaContext,
    /** The origin calls through this object run as. */
    readonly origin: RunOrigin,
    /** The group calls through this object belong to (`api.group`). */
    protected readonly group?: object,
  ) {}

  protected operation<P extends PD.Params, R>(def: OperationDefinition<P, R>): OperationFn<P, R> {
    this.ctx.runner.define(def as OperationDefinition<any, any>, this.constructor);
    const fn = ((values?: Partial<PD.ValuesFor<P>>, options?: RunOptions) =>
      this.ctx.runner.run(def, values ?? {}, {
        origin: this.origin,
        ...(this.group ? { group: this.group } : {}),
        ...options,
      })) as OperationFn<P, R>;
    Object.defineProperty(fn, 'definition', { value: def, enumerable: true });
    return fn;
  }

  /** The operations this area declares. */
  operations(): OperationFn<any, any>[] {
    return Object.values(this).filter((v): v is OperationFn<any, any> => typeof v === 'function' && 'definition' in v);
  }
}

/** The core areas (api.md §2), filled in as each is built. */
export interface CoreApi {
  readonly canvas: import('./canvas.js').CanvasApi;
  readonly camera: import('./camera.js').CameraApi;
  readonly load: import('./load.js').LoadApi;
  readonly session: import('./session.js').SessionApi;
  readonly structure: import('./structure.js').StructureApi;
  readonly history: import('./history.js').HistoryApi;
  readonly view: import('./view.js').ViewApi;
}

/** Extensions augment this with their areas, keyed by their id (api.md §6). */
export interface ExtensionApis {}

/** Every operation, as commands, tools, scripts and view models call them. */
export type ViewerApi = CoreApi &
  Partial<ExtensionApis> & {
    /** Runs an operation by id (commands, tools and the run log name operations by id). */
    run(id: string, values?: Readonly<Record<string, unknown>>, options?: RunOptions): Promise<ApiResult>;
    /** The same API bound to another origin. */
    as(origin: RunOrigin): ViewerApi;
    /**
     * The calls inside are one history step and one record group; a failure (a throw, or a failed result returned)
     * restores the group's checkpoint: all or nothing (api.md §5).
     */
    group<T>(label: string, work: (api: ViewerApi) => Promise<T>): Promise<ApiResult<T>>;
    /** An area by id, core or an extension's. */
    area(id: string): ApiArea | undefined;
    /** The runs in progress. */
    readonly running: ReadonlyState<readonly RunRecord[]>;
    readonly events: {
      readonly started: Observable<RunRecord>;
      readonly finished: Observable<RunRecord>;
    };
  };
