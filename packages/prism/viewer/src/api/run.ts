/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import { observableState, type OptimisticProjection, type ReadonlyState } from '@molstar/prism-ui';
import { type Observable, Subject } from 'rxjs';
import {
  type ApiContext,
  type ApiFailure,
  type ApiModules,
  ApiResult,
  type OperationDefinition,
  type RunOrigin,
  type ViewerApi,
} from './area.js';
import type { ViewerManagers } from './managers.js';
import type { TargetResolver } from './targets.js';
import type { PluginContext } from '@molstar/plugin/context';
import type { ViewerSources } from '../sources/viewer-sources.js';

// The one run path (operations.md §1, api.md §5): every call, from a command, a tool, a script, a URL or the
// sequence, goes through `ApiRunner.run`. It fills the params' defaults and checks the values, checks the host's modules
// and `when`, says "already so" for what exists (decision 30), opens a history step for an undoable operation (the
// checkpoint before it; a failure restores it), applies the projection, runs, and records the run.

/** One run, as the record and History read it. */
export interface RunRecord {
  readonly id: number;
  readonly operation: string;
  readonly title: string;
  readonly values: Readonly<Record<string, unknown>>;
  readonly origin: RunOrigin;
  readonly started: number;
  /** Set when it finished. */
  readonly result?: ApiResult;
  readonly ms?: number;
  /** The history step it is part of (its id), when it made one. */
  readonly step?: number;
  /** The group it ran in (`api.group`, `viewer.run`), by label. */
  readonly group?: string;
}

export interface RunInput {
  readonly origin: RunOrigin;
  readonly signal?: AbortSignal;
  readonly copy?: boolean;
  /** The atomic group the call belongs to (an `api.group`'s bound API). */
  readonly group?: object;
}

export interface ApiRunnerDeps {
  readonly plugin: PluginContext;
  readonly sources: ViewerSources;
  readonly targets: TargetResolver;
  readonly modules: ApiModules;
  readonly managers: ViewerManagers;
  /** The API bound to an origin (and group), for `ctx.api`. */
  readonly bind: (origin: RunOrigin, group?: object) => ViewerApi;
}

export class ApiRunner {
  readonly #definitions = new Map<string, OperationDefinition<any, any>>();
  readonly #owners = new Map<string, unknown>();
  readonly #projections = new Map<string, OptimisticProjection<unknown>>();
  readonly #running = observableState<readonly RunRecord[]>([]);
  readonly #started = new Subject<RunRecord>();
  readonly #finished = new Subject<RunRecord>();
  readonly #groups = new WeakMap<object, string>();
  #deps: ApiRunnerDeps | undefined;
  #serial = 0;

  /** Late: the API is made with the runner, and the runner needs it to bind `ctx.api`. */
  connect(deps: ApiRunnerDeps): void {
    this.#deps = deps;
  }

  get running(): ReadonlyState<readonly RunRecord[]> {
    return this.#running.asReadonly;
  }

  get started(): Observable<RunRecord> {
    return this.#started.asObservable();
  }

  get finished(): Observable<RunRecord> {
    return this.#finished.asObservable();
  }

  /** Every operation declared (by any area): for `api.run(id)`, commands, tools and the coverage check. */
  get definitions(): ReadonlyMap<string, OperationDefinition<any, any>> {
    return this.#definitions;
  }

  /**
   * Declares an operation. An area declares its operations again for each origin it is bound to (the same owner):
   * the first stays. Another owner under the same id throws.
   */
  define(def: OperationDefinition<any, any>, owner: unknown): void {
    const existing = this.#owners.get(def.id);
    if (existing !== undefined) {
      if (existing !== owner) throw new Error(`Operation '${def.id}' is defined twice`);
      return;
    }
    this.#owners.set(def.id, owner);
    this.#definitions.set(def.id, def);
  }

  /** Removes an operation (its extension was removed). */
  undefine(id: string): void {
    this.#definitions.delete(id);
    this.#owners.delete(id);
  }

  /** A view model's optimistic projection, which operations' `project` contribute to by its id. */
  registerProjection(model: string, projection: OptimisticProjection<unknown>): () => void {
    this.#projections.set(model, projection);
    return () => {
      if (this.#projections.get(model) === projection) this.#projections.delete(model);
    };
  }

  async run<P extends PD.Params, R>(
    def: OperationDefinition<P, R>,
    input: Partial<PD.ValuesFor<P>> | Readonly<Record<string, unknown>>,
    run: RunInput,
  ): Promise<ApiResult<R>> {
    const deps = this.#deps;
    if (!deps) return ApiResult.fail('unavailable', 'The viewer is not ready yet.') as ApiResult<R>;
    const signal = run.signal ?? new AbortController().signal;
    const ctx: ApiContext = {
      plugin: deps.plugin,
      sources: deps.sources,
      targets: deps.targets,
      api: deps.bind(run.origin, run.group),
      origin: run.origin,
      modules: deps.modules,
      managers: deps.managers,
      signal,
      copy: !!run.copy,
    };

    const params = typeof def.params === 'function' ? def.params(ctx) : def.params;
    const checked = checkValues(params, input as Readonly<Record<string, unknown>>);
    if ('ok' in checked) return this.#refused(def, input, run, checked) as ApiResult<R>;
    const values = checked.values as PD.ValuesFor<P>;

    for (const m of def.modules ?? []) {
      if (!deps.modules.has(m)) {
        return this.#refused(
          def,
          values,
          run,
          ApiResult.fail('unavailable', `${def.title} needs ${ModuleNames[m]}, which this viewer does not have.`),
        ) as ApiResult<R>;
      }
    }
    const when = def.when?.(ctx) ?? true;
    if (when !== true) return this.#refused(def, values, run, ApiResult.fail('blocked', when)) as ApiResult<R>;
    if (!ctx.copy) {
      const exists = def.exists?.(values, ctx);
      if (exists) return this.#refused(def, values, run, ApiResult.ok(exists, {})) as ApiResult<R>;
    }
    if (signal.aborted) return this.#refused(def, values, run, aborted()) as ApiResult<R>;

    const record: RunRecord = {
      id: ++this.#serial,
      operation: def.id,
      title: def.title,
      values: values as Record<string, unknown>,
      origin: run.origin,
      started: Date.now(),
      ...(run.group && this.#groups.has(run.group) ? { group: this.#groups.get(run.group) } : {}),
    };
    this.#running.set([...this.#running.value, record]);
    this.#started.next(record);

    const projection = def.project ? this.#project(def.project(values)) : undefined;
    const work = async (): Promise<ApiResult<R>> => {
      try {
        if (signal.aborted) return aborted() as ApiResult<R>;
        const result = await def.run(values, ctx);
        return signal.aborted && !result.ok ? (aborted() as ApiResult<R>) : result;
      } catch (e) {
        if (signal.aborted || isAbort(e)) return aborted() as ApiResult<R>;
        return ApiResult.fail('failed', `${def.title} failed: ${messageOf(e)}`) as ApiResult<R>;
      }
    };

    const t0 = performance.now();
    let result: ApiResult<R>;
    let step: number | undefined;
    if (def.read || def.history === false) {
      result = await work();
    } else {
      try {
        const outcome = await deps.managers.history.run(
          (r: ApiResult<R>) => (r.ok ? stepLabel(def.title, r.did) : def.title),
          run.origin,
          work,
          { failed: (r) => !r.ok, ...(def.undoable === false ? { undoable: false } : {}) },
        );
        result = outcome.value;
        step = outcome.step?.id;
        if (!result.ok && outcome.restored) {
          result = { ...result, hint: result.hint ?? 'Nothing was changed: the viewer is as it was before.' };
        }
      } catch (e) {
        result = ApiResult.fail('failed', `${def.title} failed: ${messageOf(e)}`) as ApiResult<R>;
      }
    }
    if (projection) {
      if (result.ok) projection.acknowledge();
      else projection.reject();
    }
    const done: RunRecord = {
      ...record,
      result,
      ms: Math.round(performance.now() - t0),
      ...(step !== undefined ? { step } : {}),
    };
    this.#running.set(this.#running.value.filter((r) => r.id !== record.id));
    this.#finished.next(done);
    return result;
  }

  /**
   * Raw plugin work as one history step (`viewer.run`, README §3.2): recorded as one opaque run, `script.run`, with
   * its label. A throw restores the checkpoint and is thrown on.
   */
  async opaque<T>(label: string, origin: RunOrigin, work: () => Promise<T>): Promise<T> {
    const deps = this.#deps;
    if (!deps) throw new Error('The viewer is not ready yet.');
    const record: RunRecord = {
      id: ++this.#serial,
      operation: 'script.run',
      title: label,
      values: {},
      origin,
      started: Date.now(),
    };
    this.#running.set([...this.#running.value, record]);
    this.#started.next(record);
    const t0 = performance.now();
    let result: ApiResult = ApiResult.ok(label);
    let step: number | undefined;
    try {
      const outcome = await deps.managers.history.run(label, origin, work);
      step = outcome.step?.id;
      return outcome.value;
    } catch (e) {
      result = ApiResult.fail('failed', `${label} failed: ${messageOf(e)}`, 'Nothing was changed.');
      throw e;
    } finally {
      this.#running.set(this.#running.value.filter((r) => r.id !== record.id));
      this.#finished.next({
        ...record,
        result,
        ms: Math.round(performance.now() - t0),
        ...(step !== undefined ? { step } : {}),
      });
    }
  }

  /** Names an atomic group's token, for the record. */
  nameGroup(token: object, label: string): void {
    this.#groups.set(token, label);
  }

  dispose(): void {
    this.#started.complete();
    this.#finished.complete();
    this.#definitions.clear();
    this.#projections.clear();
  }

  #project(p: { model: string; apply(base: unknown): unknown; isAbsorbed(base: unknown): boolean }) {
    const projection = this.#projections.get(p.model);
    const deps = this.#deps;
    if (!projection || !deps) return undefined;
    const tx = deps.managers.transactions.begin();
    projection.contribute(tx, { apply: p.apply, isAbsorbed: p.isAbsorbed });
    return tx;
  }

  /** A call that did not run: recorded too, so the record says what was asked and why nothing happened. */
  #refused(
    def: OperationDefinition<any, any>,
    values: Readonly<Record<string, unknown>>,
    run: RunInput,
    result: ApiResult,
  ): ApiResult {
    const record: RunRecord = {
      id: ++this.#serial,
      operation: def.id,
      title: def.title,
      values,
      origin: run.origin,
      started: Date.now(),
      result,
      ms: 0,
    };
    this.#finished.next(record);
    return result;
  }
}

const ModuleNames: Record<string, string> = {
  canvas: 'a canvas',
  files: 'access to files',
  network: 'the network',
  storage: 'storage',
  page: 'a page',
};

function aborted(): ApiFailure {
  return ApiResult.fail('aborted', 'Stopped before it finished.', 'Nothing was changed.');
}

function isAbort(e: unknown): boolean {
  return (
    !!e &&
    typeof e === 'object' &&
    ((e as { name?: string }).name === 'AbortError' || (e as { isAborted?: boolean }).isAborted === true)
  );
}

function messageOf(e: unknown): string {
  if (e instanceof Error) return e.message;
  if (typeof e === 'string') return e;
  try {
    return JSON.stringify(e);
  } catch {
    return String(e);
  }
}

/** A history step's label: what the result says it did, without its full stop; the title when it says nothing. */
function stepLabel(title: string, did: string): string {
  const text = did.trim().replace(/\.$/, '');
  return text || title;
}

/**
 * Fills a call's values from the params' defaults and checks them: unknown names, wrong types, a select's options, a
 * number's range. Nested groups are filled the same way. A failure names the field and what it takes.
 */
export function checkValues(
  params: PD.Params,
  input: Readonly<Record<string, unknown>>,
): { values: Record<string, unknown> } | ApiFailure {
  const values: Record<string, unknown> = PD.getDefaultValues(params) as Record<string, unknown>;
  for (const [name, value] of Object.entries(input ?? {})) {
    if (value === undefined) continue;
    const p = params[name];
    if (!p) {
      const known = Object.keys(params);
      return ApiResult.fail(
        'invalid',
        `Unknown value '${name}'.`,
        known.length ? `It takes: ${known.join(', ')}.` : 'It takes no values.',
      );
    }
    const problem = checkParam(name, p, value);
    if (problem) return problem;
    values[name] =
      p.type === 'group' && value && typeof value === 'object'
        ? { ...(values[name] as object), ...(value as object) }
        : value;
  }
  return { values };
}

function checkParam(name: string, p: PD.Any, value: unknown): ApiFailure | undefined {
  switch (p.type) {
    case 'number':
    case 'converted': {
      if (p.type === 'converted') return undefined;
      if (typeof value !== 'number' || Number.isNaN(value))
        return ApiResult.fail('invalid', `'${name}' takes a number.`);
      if ((p.min !== undefined && value < p.min) || (p.max !== undefined && value > p.max)) {
        return ApiResult.fail('invalid', `'${name}' takes a number from ${p.min ?? '…'} to ${p.max ?? '…'}.`);
      }
      return undefined;
    }
    case 'boolean':
      return typeof value === 'boolean' ? undefined : ApiResult.fail('invalid', `'${name}' takes true or false.`);
    case 'text':
      return typeof value === 'string' ? undefined : ApiResult.fail('invalid', `'${name}' takes text.`);
    case 'select': {
      const options = p.options.map((o) => o[0]);
      return options.includes(value as never)
        ? undefined
        : ApiResult.fail(
            'invalid',
            `'${name}' cannot be ${JSON.stringify(value)}.`,
            `It takes: ${options.join(', ')}.`,
          );
    }
    case 'group':
      return value && typeof value === 'object' && !Array.isArray(value)
        ? undefined
        : ApiResult.fail('invalid', `'${name}' takes an object.`);
    default:
      return undefined;
  }
}
