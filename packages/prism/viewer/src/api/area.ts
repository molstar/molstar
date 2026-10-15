/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import type { PluginContext } from '@molstar/plugin/context';
import type { HistoryContext, ReadonlyState } from '@molstar/prism-ui/core';
import type { Observable } from 'rxjs';
import type { ViewerSources } from '../sources/viewer-sources.js';
import type { ViewerManagers } from './managers.js';
import type { RunRecord } from './run.js';
import type { TargetResolver } from './targets.js';
import type { ExtensionApis } from '../api.js';

// The viewer's API (api.md, operations.md §1): one object per area, its operations as fields, every change through
// one run path (`ApiRunner`), results for every caller, never a throw. An area is made once; binding it to an origin
// (`api.as`), a group (`api.group`) or a run in progress (`ctx.api`) is a thin facade over the same object.

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

/**
 * Why a call did not do what was asked: `invalid` (its values), `blocked` (not now: `when`), `unavailable` (the host
 * lacks a module), `not-found`, `ambiguous` (which one?), `confirm` (ask first, then call again), `aborted`, `failed`.
 */
export type ApiErrorCode =
  | 'invalid'
  | 'blocked'
  | 'unavailable'
  | 'not-found'
  | 'ambiguous'
  | 'confirm'
  | 'aborted'
  | 'failed';

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
      /** It was already so: nothing changed, no history step (operations.md §1, desired state). */
      readonly unchanged?: true;
    }
  | {
      readonly ok: false;
      readonly code: ApiErrorCode;
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
  /** It is already so: nothing changed (no history step). */
  already<R>(did: string, extra: { value?: R } = {}): ApiResult<R> {
    return { ok: true, did, changed: {}, unchanged: true, ...extra };
  },
  fail(code: ApiErrorCode, message: string, hint?: string): ApiFailure {
    return hint ? { ok: false, code, message, hint } : { ok: false, code, message };
  },
  isFailure(x: unknown): x is ApiFailure {
    return !!x && typeof x === 'object' && (x as { ok?: unknown }).ok === false && 'code' in (x as object);
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

/** What an area and an operation's `params` and `when` get: the viewer's parts, no origin (api.md §4). */
export interface ApiAreaContext {
  /** Only `api/` and `sources/` use it (contexts.md §7). */
  readonly plugin: PluginContext;
  readonly sources: ViewerSources;
  readonly targets: TargetResolver;
  readonly modules: ApiModules;
  readonly managers: ViewerManagers;
  /** Declares and runs operations: the run path, as an area sees it. */
  readonly operations: OperationPort;
}

/** What an operation's run gets: the area context, bound to who runs it (api.md §4). */
export interface ApiContext extends ApiAreaContext {
  /** Other areas, bound to the same origin, and part of this run: calls through it join its history step. */
  readonly api: ViewerApi;
  readonly origin: RunOrigin;
  readonly signal: AbortSignal;
  /** Asked for another of what exists (decision 30). */
  readonly copy: boolean;
}

/** The run path as an area sees it. */
export interface OperationPort {
  /** Declares an operation of an area; its id must start with the area's id (`<area>.<verb>`). */
  define(def: OperationDefinition<any, any>, area: string): void;
  run<P extends PD.Params, R>(
    def: OperationDefinition<P, R>,
    values: Partial<PD.ValuesFor<P>> | Readonly<Record<string, unknown>>,
    input: RunInput,
  ): Promise<ApiResult<R>>;
}

/** How a call runs: who, and what it is part of. */
export interface RunInput {
  readonly origin: RunOrigin;
  readonly signal?: AbortSignal;
  readonly copy?: boolean;
  /** The run or group in progress the call is part of (no step of its own). */
  readonly within?: HistoryContext;
  /** The group's record, which collects its calls. */
  readonly group?: GroupRecord;
}

/** A group's record while it runs: its id and the runs made inside it. */
export interface GroupRecord {
  readonly id: number;
  readonly label: string;
  readonly runs: RunRecord[];
}

export interface OperationDefinition<P extends PD.Params = PD.Params, R = unknown> {
  /** `<area>.<verb>`: `component.add`, `canvas.background`. */
  readonly id: string;
  /** For people: the command's title, and the history entry when the result says nothing better. */
  readonly title: string;
  /** One definition: the command's form and the tool's schema derive from it; a function for registry-driven ones. */
  readonly params: P | ((ctx: ApiAreaContext) => P);
  /** Changes nothing: no history step, may run alongside others. */
  readonly read?: true;
  /** Changes view state, not content: no history step (select, focus, camera, highlight, view). */
  readonly history?: false;
  /** Changes nothing the history could restore (an upload): no checkpoint, no step; default true. */
  readonly undoable?: false;
  /** Asks before it runs, from the agent (agent.md §4). */
  readonly approval?: 'network' | 'upload' | 'destructive';
  /** What it needs of the host; where one is missing it says why instead of running. */
  readonly modules?: readonly ApiModuleKey[];
  /** Whether it can run now; the string is the reason. Commands and tools read it too, so no origin. */
  when?(ctx: ApiAreaContext): true | string;
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

export function isOperationFn(x: unknown): x is OperationFn<any, any> {
  return typeof x === 'function' && 'definition' in x;
}

const Script: RunOrigin = { kind: 'script' };

/**
 * The base of every area: its id, and `operation` to declare an operation as a field (api.md §3). Made once per viewer;
 * its operation fields run as `script`, and the API's facades rebind them to other origins.
 */
export abstract class ApiArea {
  abstract readonly id: string;

  constructor(protected readonly ctx: ApiAreaContext) {}

  protected operation<P extends PD.Params, R>(def: OperationDefinition<P, R>): OperationFn<P, R> {
    // `id` is declared before the operation fields, so it is set by now.
    this.ctx.operations.define(def as OperationDefinition<any, any>, this.id);
    return bindOperation(this.ctx.operations, def, { origin: Script });
  }

  /**
   * Adds an operation to this area from outside it (a built-in feature's: `loaders` brings `load.entry`), under
   * `name`; returns its removal. The id must start with this area's id.
   */
  contribute<P extends PD.Params, R>(name: string, def: OperationDefinition<P, R>): () => void {
    if (name in this) throw new Error(`Area '${this.id}' already has '${name}'`);
    this.ctx.operations.define(def as OperationDefinition<any, any>, this.id);
    const fn = bindOperation(this.ctx.operations, def, { origin: Script });
    Object.defineProperty(this, name, { value: fn, enumerable: true, configurable: true });
    return () => {
      const self = this as unknown as Record<string, unknown>;
      if (self[name] === fn) delete self[name];
    };
  }

  /** The operations this area declares, contributed ones included. */
  operations(): OperationFn<any, any>[] {
    return Object.values(this).filter(isOperationFn);
  }

  /** Lets go of what the area holds (an extension's area when the extension goes). */
  dispose?(): void;
}

/** An operation bound to how it runs. */
export function bindOperation<P extends PD.Params, R>(
  port: OperationPort,
  def: OperationDefinition<P, R>,
  input: RunInput,
): OperationFn<P, R> {
  const fn = ((values?: Partial<PD.ValuesFor<P>>, options?: RunOptions) =>
    port.run(def, values ?? {}, { ...input, ...options })) as OperationFn<P, R>;
  Object.defineProperty(fn, 'definition', { value: def, enumerable: true });
  return fn;
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

/** An operation as commands and tools see it, before running: its params resolved, and whether it can run now. */
export interface OperationDescription {
  readonly definition: OperationDefinition<any, any>;
  readonly params: PD.Params;
  /** True, or why not. */
  readonly available: true | string;
}

export interface GroupOptions {
  readonly signal?: AbortSignal;
  /**
   * A call inside that fails: `rollback` (default) fails the whole group at its end, restoring its checkpoint;
   * `continue` undoes just that call and goes on.
   */
  readonly onFailure?: 'rollback' | 'continue';
}

/** A group's result: the work's value, what changed (the union of its calls'), and its calls. */
export type GroupResult<T> = ApiResult<T> & { readonly runs: readonly RunRecord[] };

/** Every operation, as commands, tools, scripts and view models call them. */
export type ViewerApi = CoreApi &
  Partial<ExtensionApis> & {
    /** Runs an operation by id (commands, tools and the run log name operations by id). */
    run(id: string, values?: Readonly<Record<string, unknown>>, options?: RunOptions): Promise<ApiResult>;
    /** The same API bound to another origin (inside a group or a run, still part of it). */
    as(origin: RunOrigin): ViewerApi;
    /**
     * The calls inside are one history step and one record group; a failure (a throw, a failed result returned, or a
     * call inside that failed) restores the group's checkpoint: all or nothing (api.md §5).
     */
    group<T>(label: string, work: (api: ViewerApi) => Promise<T>, options?: GroupOptions): Promise<GroupResult<T>>;
    /** An area by id, core or an extension's (bound like this API). */
    area(id: string): ApiArea | undefined;
    /** An operation by id, as commands and tools see it. */
    describe(id: string): OperationDescription | undefined;
    /** The operations declared, by id; changes when an extension adds or removes some. */
    readonly definitions: ReadonlyState<ReadonlyMap<string, OperationDefinition<any, any>>>;
    /** The runs in progress. */
    readonly running: ReadonlyState<readonly RunRecord[]>;
    readonly events: {
      readonly started: Observable<RunRecord>;
      readonly finished: Observable<RunRecord>;
    };
  };
