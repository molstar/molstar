/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import { type Observable } from 'rxjs';
import { type ObservableState, observableState } from '../reactive/observable-state.js';
import type { ReadonlyState } from '../reactive/readonly-state.js';

import { derive } from '../reactive/state.js';
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
  readonly #folds = new Map<string, ObservableState<boolean>>();

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
      stream = derive(this.#values.asReadonly, (v) => v[name]);
      this.#fields.set(name, stream);
    }
    return stream as Observable<PD.Values<P>[K]>;
  }

  /** Is the field shown now (its `hideIf` evaluated on the current values)? */
  visible(name: keyof P & string): Observable<boolean> {
    const param = this.params[name];
    if (!param.hideIf) return derive(this.#values.asReadonly, () => true);
    return derive(this.#values.asReadonly, (v) => !param.hideIf!(v));
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

  /**
   * Whether a fold of the form is open (Advanced, a category, a choice's params), by its path in the form: kept for the
   * form's life, so a fold stays as the user left it when its rows render again (another choice and back, a reload).
   */
  fold(path: string, initial = false): ObservableState<boolean> {
    let state = this.#folds.get(path);
    if (!state) {
      state = observableState(initial);
      this.#folds.set(path, state);
    }
    return state;
  }
}

/** Params under one heading: a category's (PD `category`; its own fold), or those without one. */
export interface ParamSection {
  category?: string;
  names: string[];
}

/**
 * The params as the classic UI lays them out (v5 `ParameterControls`): the essential ones (`isEssential`) and the rest,
 * each in sections by category (those without one first, then the categories by name); hidden params and plain values
 * (`PD.Value`, data carried along) have no row. The rest fold under "Advanced" only when there are essential ones too.
 */
export function paramSections(params: PD.Params): { essentials: ParamSection[]; advanced: ParamSection[] } {
  const split = (essential: boolean) => {
    const plain: ParamSection = { names: [] };
    const byCategory = new Map<string, ParamSection>();
    for (const name of Object.keys(params)) {
      const param = params[name];
      if (param.isHidden || param.type === 'value' || !!param.isEssential !== essential) continue;
      if (!param.category) {
        plain.names.push(name);
        continue;
      }
      let section = byCategory.get(param.category);
      if (!section) {
        section = { category: param.category, names: [] };
        byCategory.set(param.category, section);
      }
      section.names.push(name);
    }
    const categories = [...byCategory.values()].sort((a, b) => (a.category! < b.category! ? -1 : 1));
    return plain.names.length ? [plain, ...categories] : categories;
  };
  return { essentials: split(true), advanced: split(false) };
}
