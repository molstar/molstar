/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { isBytes, takesBytes } from './input.js';
import { isNodeQuery } from './nodes.js';
import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import type { PluginContext } from '@molstar/plugin/context';
import {
  type HistoryContext,
  observableState,
  type OptimisticProjection,
  type ReadonlyState,
} from '@molstar/prism-ui/core';
import { type Observable, Subject } from 'rxjs';
import type { ViewerSources } from '../sources/viewer-sources.js';
import {
  type ApiAreaContext,
  type GroupRecord,
  type ApiContext,
  type ApiFailure,
  type ApiModules,
  ApiResult,
  type OperationDefinition,
  type OperationDescription,
  type OperationPort,
  type RunInput,
  type RunOrigin,
  type ViewerApi,
} from './area.js';
import type { ViewerManagers } from './managers.js';
import type { TargetResolver } from './targets.js';

// The one run path (operations.md §1, api.md §5): every call, from a command, a tool, a script, a URL or the
// sequence, goes through `ApiRunner.run`. It fills the params' defaults and checks the values, checks the host's modules
// and `when`, says "already so" for what exists (decision 30), runs an undoable operation as a history step (the
// checkpoint before it; a failure restores it; a result that changed nothing records none), applies the projection,
// and records the run. A call made inside a run or a group is part of it: its `ctx.api` carries the history context.
// It never throws: a failure is a result.

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
  /** The history step it made or joined (its id). */
  readonly step?: number;
  /** The group it ran in (`api.group`), by id and label. */
  readonly group?: { readonly id: number; readonly label: string };
  /** It was part of another run (an operation calling another). */
  readonly nested?: true;
}

export interface ApiRunnerDeps {
  readonly plugin: PluginContext;
  readonly sources: ViewerSources;
  readonly targets: TargetResolver;
  readonly modules: ApiModules;
  readonly managers: ViewerManagers;
  /** The API bound to how a call runs, for `ctx.api`. */
  readonly bind: (input: RunInput) => ViewerApi;
}

export class ApiRunner implements OperationPort {
  readonly #aliases = new Map<string, string>();
  readonly #definitions = observableState<ReadonlyMap<string, OperationDefinition<any, any>>>(new Map(), {
    compare: Object.is,
  });
  readonly #owners = new Map<string, string>();
  readonly #projections = new Map<string, OptimisticProjection<unknown>>();
  readonly #running = observableState<readonly RunRecord[]>([]);
  readonly #started = new Subject<RunRecord>();
  readonly #finished = new Subject<RunRecord>();
  #deps: ApiRunnerDeps | undefined;
  #areaContext: ApiAreaContext | undefined;
  #serial = 0;
  #disposed = false;

  /** Late: the API is made with the runner, and the runner needs it to bind `ctx.api`. */
  connect(deps: ApiRunnerDeps, areaContext: ApiAreaContext): void {
    this.#deps = deps;
    this.#areaContext = areaContext;
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

  /** Every operation declared, by id: for `api.run(id)`, commands, tools and the coverage check. */
  get definitions(): ReadonlyState<ReadonlyMap<string, OperationDefinition<any, any>>> {
    return this.#definitions.asReadonly;
  }

  /**
   * Declares an operation of an area. Its id must start with the area's id. An area declares its operations once;
   * another area declaring an id that exists throws.
   */
  define(def: OperationDefinition<any, any>, area: string): void {
    if (!def.id.startsWith(`${area}.`)) {
      throw new Error(`Operation '${def.id}' must start with its area's id, '${area}.'`);
    }
    if (this.#owners.has(def.id)) {
      if (this.#definitions.value.get(def.id) === def) return;
      throw new Error(`Operation '${def.id}' is defined twice`);
    }
    this.#owners.set(def.id, area);
    // An alias may name another area (`structure.model` runs `trajectory.frame`): one id, one operation.
    for (const alias of def.aliases ?? []) {
      if (this.#owners.has(alias) || this.#aliases.has(alias)) throw new Error(`Alias '${alias}' is taken`);
      this.#aliases.set(alias, def.id);
    }
    this.#definitions.set(new Map([...this.#definitions.value, [def.id, def]]));
  }

  /** An operation by id or by a former id (`aliases`). */
  definition(id: string): OperationDefinition<any, any> | undefined {
    const defs = this.#definitions.value;
    return defs.get(id) ?? defs.get(this.#aliases.get(id) ?? '');
  }

  /** Removes an operation (its extension was removed, or a stub was replaced). */
  undefine(id: string): void {
    if (!this.#owners.delete(id)) return;
    for (const [alias, to] of this.#aliases) if (to === id) this.#aliases.delete(alias);
    const next = new Map(this.#definitions.value);
    next.delete(id);
    this.#definitions.set(next);
  }

  /** An operation as commands and tools see it: params resolved, whether it can run now. */
  describe(id: string): OperationDescription | undefined {
    const definition = this.definition(id);
    const ctx = this.#areaContext;
    if (!definition || !ctx) return undefined;
    try {
      const params = typeof definition.params === 'function' ? definition.params(ctx) : definition.params;
      return { definition, params, available: availability(definition, ctx), supported: support(definition, ctx) };
    } catch (e) {
      const why = `${definition.title} cannot run: ${messageOf(e)}`;
      return { definition, params: {}, available: why, supported: why };
    }
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
    const area = this.#areaContext;
    if (!deps || !area || this.#disposed) {
      return ApiResult.fail(
        'unavailable',
        this.#disposed ? 'The viewer is gone.' : 'The viewer is not ready yet.',
      ) as ApiResult<R>;
    }
    const signal = run.signal ?? new AbortController().signal;
    const history = deps.managers.history;
    // A context from a run that has ended is stale: the call runs on its own.
    const within = history.isActive(run.within) ? run.within : undefined;
    const makeCtx = (context: HistoryContext | undefined): ApiContext => ({
      ...area,
      api: deps.bind({ ...run, ...(context ? { within: context } : {}) }),
      origin: run.origin,
      signal,
      copy: !!run.copy,
    });

    let values: PD.ValuesFor<P>;
    try {
      const params = typeof def.params === 'function' ? def.params(area) : def.params;
      const checked = def.passValues
        ? { values: { ...(input as Record<string, unknown>) } }
        : checkValues(params, input as Readonly<Record<string, unknown>>);
      if ('ok' in checked) return this.#refused(def, input, run, checked) as ApiResult<R>;
      values = checked.values as PD.ValuesFor<P>;
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
      const when = def.when?.(area) ?? true;
      if (when !== true) return this.#refused(def, values, run, ApiResult.fail('blocked', when)) as ApiResult<R>;
      if (!run.copy) {
        const exists = def.exists?.(values, makeCtx(within));
        if (exists) return this.#refused(def, values, run, ApiResult.already(exists)) as ApiResult<R>;
      }
    } catch (e) {
      return this.#refused(
        def,
        input,
        run,
        ApiResult.fail('failed', `${def.title} failed: ${messageOf(e)}`),
      ) as ApiResult<R>;
    }
    if (signal.aborted) return this.#refused(def, values, run, aborted()) as ApiResult<R>;

    const record: RunRecord = {
      id: ++this.#serial,
      operation: def.id,
      title: def.title,
      values: values as Record<string, unknown>,
      origin: run.origin,
      started: Date.now(),
      ...(run.group ? { group: { id: run.group.id, label: run.group.label } } : {}),
      ...(within ? { nested: true as const } : {}),
    };
    if (run.group && within && within === run.group.within) run.group.direct.add(record.id);
    this.#running.set([...this.#running.value, record]);
    this.#started.next(record);

    const projection = def.project ? this.#project(def.project(values)) : undefined;
    const work = async (context: HistoryContext | undefined): Promise<ApiResult<R>> => {
      try {
        if (signal.aborted) return aborted() as ApiResult<R>;
        const result = await def.run(values, makeCtx(context));
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
      result = await work(within);
    } else {
      try {
        const outcome = await history.run(
          (r: ApiResult<R>) => (r.ok ? stepLabel(def.title, r.did) : def.title),
          run.origin,
          work,
          {
            failed: (r) => !r.ok,
            record: (r) => r.ok && !r.unchanged,
            ...(def.undoable === false ? { undoable: false } : {}),
            ...(within ? { within } : {}),
          },
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
    run.group?.runs.push(done);
    if (!this.#disposed) this.#finished.next(done);
    return result;
  }

  /**
   * Raw plugin work as one history step (`viewer.run`, README §3.2): recorded as one opaque run, `script.run`, with
   * its label. Calls it makes through the API with its context join the step. A throw restores the checkpoint and is
   * thrown on.
   */
  async opaque<T>(label: string, origin: RunOrigin, work: (context: HistoryContext) => Promise<T>): Promise<T> {
    const deps = this.#deps;
    if (!deps || this.#disposed) throw new Error('The viewer is not ready, or is gone.');
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
      if (!this.#disposed) {
        this.#finished.next({
          ...record,
          result,
          ms: Math.round(performance.now() - t0),
          ...(step !== undefined ? { step } : {}),
        });
      }
    }
  }

  /** A new group's record. */
  newGroup(label: string): GroupRecord {
    return { id: ++this.#serial, label, runs: [], direct: new Set() };
  }

  /** No more runs: calls after this say the viewer is gone. */
  dispose(): void {
    this.#disposed = true;
    this.#started.complete();
    this.#finished.complete();
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
      ...(run.group ? { group: { id: run.group.id, label: run.group.label } } : {}),
    };
    run.group?.runs.push(record);
    if (!this.#disposed) this.#finished.next(record);
    return result;
  }
}

function support(def: OperationDefinition<any, any>, ctx: ApiAreaContext): true | string {
  for (const m of def.modules ?? []) {
    if (!ctx.modules.has(m)) return `${def.title} needs ${ModuleNames[m]}, which this viewer does not have.`;
  }
  return true;
}

function availability(def: OperationDefinition<any, any>, ctx: ApiAreaContext): true | string {
  const supported = support(def, ctx);
  return supported === true ? (def.when?.(ctx) ?? true) : supported;
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

export function messageOf(e: unknown): string {
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
    // A node(...) query is its selector text (nodes.ts).
    if (p?.type === 'text' && isNodeQuery(value)) {
      values[name] = value.text;
      continue;
    }
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
    case 'number': {
      if (typeof value !== 'number' || Number.isNaN(value)) {
        return ApiResult.fail('invalid', `'${name}' takes a number.`);
      }
      if ((p.min !== undefined && value < p.min) || (p.max !== undefined && value > p.max)) {
        return ApiResult.fail('invalid', `'${name}' takes a number from ${p.min ?? '…'} to ${p.max ?? '…'}.`);
      }
      return undefined;
    }
    case 'boolean':
      return typeof value === 'boolean' ? undefined : ApiResult.fail('invalid', `'${name}' takes true or false.`);
    case 'text':
      // `data` (api/input.ts) takes bytes as well.
      if (typeof value === 'string' || (takesBytes(p) && isBytes(value))) return undefined;
      return ApiResult.fail('invalid', `'${name}' takes text${takesBytes(p) ? ' or bytes' : ''}.`);
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
    case 'mapped': {
      // A switch or a choice with settings (operations.md §1): `{ name, params }`, the params the named case's.
      const names = p.select.options.map((o) => o[0]);
      const v = value as { name?: unknown; params?: unknown } | null;
      if (!v || typeof v !== 'object' || Array.isArray(v) || !names.includes(v.name as string)) {
        return ApiResult.fail(
          'invalid',
          `'${name}' takes { name, params } with a name of: ${names.join(', ')}.`,
          names.includes('on') ? `For example { name: 'on' } or { name: 'off' }.` : undefined,
        );
      }
      if (v.params === undefined) return undefined;
      const inner = p.map(v.name as string);
      if (inner.type !== 'group') return undefined;
      if (!v.params || typeof v.params !== 'object' || Array.isArray(v.params)) {
        return ApiResult.fail('invalid', `'${name}.params' takes an object.`);
      }
      for (const [k, val] of Object.entries(v.params as Record<string, unknown>)) {
        if (val === undefined) continue;
        const q = inner.params[k];
        if (!q) {
          const known = Object.keys(inner.params);
          return ApiResult.fail(
            'invalid',
            `Unknown value '${name}.params.${k}'.`,
            known.length ? `It takes: ${known.join(', ')}.` : 'It takes no values.',
          );
        }
        const problem = checkParam(`${name}.params.${k}`, q, val);
        if (problem) return problem;
      }
      return undefined;
    }
    default:
      return undefined;
  }
}
