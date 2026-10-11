/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import { distinctUntilChanged, map, type Observable } from 'rxjs';
import { observableState } from '../reactive/observable-state.js';
import type { ReadonlyState } from '../reactive/readonly-state.js';

export interface ParamsFormOptions<P extends PD.Params> {
  /** Called with the new values after each change (apply them, or keep them as a draft). */
  onChange?(values: PD.Values<P>, changed: keyof P & string): void;
}

/**
 * A form over `ParamDefinition` params: the values, a stream per field (so a control updates only when its own value
 * changes), and edits. Shared by the inspector, settings, dialogs and command forms with complex arguments.
 */
export class ParamsFormModel<P extends PD.Params = PD.Params> {
  readonly #values;
  readonly #defaults: PD.Values<P>;
  readonly #fields = new Map<string, Observable<unknown>>();

  constructor(
    readonly params: P,
    values: PD.Values<P> = PD.getDefaultValues(params),
    readonly options: ParamsFormOptions<P> = {},
  ) {
    this.#defaults = PD.getDefaultValues(params);
    this.#values = observableState<PD.Values<P>>({ ...this.#defaults, ...values }, { compare: Object.is });
  }

  get values(): ReadonlyState<PD.Values<P>> {
    return this.#values.asReadonly;
  }

  /** The value of one field, emitted when it changes (by identity). Cached per field. */
  field<K extends keyof P & string>(name: K): Observable<PD.Values<P>[K]> {
    let stream = this.#fields.get(name);
    if (!stream) {
      stream = this.#values.asReadonly.pipe(
        map((v) => v[name]),
        distinctUntilChanged(),
      );
      this.#fields.set(name, stream);
    }
    return stream as Observable<PD.Values<P>[K]>;
  }

  /** Is the field shown now (its `hideIf` evaluated on the current values)? */
  visible(name: keyof P & string): Observable<boolean> {
    const param = this.params[name];
    if (!param.hideIf)
      return this.#values.asReadonly.pipe(
        map(() => true),
        distinctUntilChanged(),
      );
    return this.#values.asReadonly.pipe(
      map((v) => !param.hideIf!(v)),
      distinctUntilChanged(),
    );
  }

  set<K extends keyof P & string>(name: K, value: PD.Values<P>[K]): void {
    const current = this.#values.value;
    if (Object.is(current[name], value)) return;
    const next = { ...current, [name]: value } as PD.Values<P>;
    this.#values.set(next);
    this.options.onChange?.(next, name);
  }

  /** Replaces all values (e.g. the authoritative ones after an update elsewhere), without reporting a change. */
  load(values: PD.Values<P>): void {
    this.#values.set({ ...this.#defaults, ...values });
  }

  /** Back to the params' defaults (reported as a change of every field that differed). */
  reset(): void {
    const current = this.#values.value;
    const next = { ...this.#defaults };
    this.#values.set(next);
    for (const name of Object.keys(this.params)) {
      if (!PD.isParamEqual(this.params[name], current[name], next[name])) {
        this.options.onChange?.(next, name as keyof P & string);
      }
    }
  }

  /** The params split as the classic UI does: essentials first; the rest under "Advanced" when there are essentials. */
  layout(): { essentials: string[]; advanced: string[] } {
    const essentials: string[] = [];
    const rest: string[] = [];
    for (const name of Object.keys(this.params)) {
      const param = this.params[name];
      if (param.isHidden) continue;
      if (param.isEssential) essentials.push(name);
      else rest.push(name);
    }
    return essentials.length ? { essentials, advanced: rest } : { essentials: rest, advanced: [] };
  }
}
