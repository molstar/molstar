/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { derive, type ReadonlyState, State } from '@molstar/core/reactive';
import { ParamDefinition as PD } from '@molstar/core/util/param-definition';

export interface ParamsFormOptions<P extends PD.Params> {
  /** Called with the new values after each change (apply them, or keep them as a draft). */
  onChange?(values: PD.Values<P>, changed: keyof P & string): void;
}

/** Shown always: a field without a `hideIf`. */
const Always: ReadonlyState<boolean> = new State(true);

/**
 * A form over `ParamDefinition` params: the values (one state), a derived state per field (so a control is told only
 * when its own value changes), and edits. Shared by inspectors, settings, dialogs and command forms with complex
 * arguments.
 */
export class ParamsFormModel<P extends PD.Params = PD.Params> {
  readonly #values: State<PD.Values<P>>;
  readonly #defaults: PD.Values<P>;
  readonly #fields = new Map<string, ReadonlyState<unknown>>();
  readonly #visible = new Map<string, ReadonlyState<boolean>>();
  readonly #folds = new Map<string, State<boolean>>();

  constructor(
    readonly params: P,
    values: PD.Values<P> = PD.getDefaultValues(params),
    readonly options: ParamsFormOptions<P> = {},
  ) {
    this.#defaults = PD.getDefaultValues(params);
    // By identity: each edit makes a new record, and an equal copy loaded from elsewhere is a change to show.
    this.#values = new State<PD.Values<P>>({ ...this.#defaults, ...values }, Object.is);
  }

  get values(): ReadonlyState<PD.Values<P>> {
    return this.#values.asReadonly;
  }

  /** The value of one field, told when it changes (by identity). Cached per field. */
  field<K extends keyof P & string>(name: K): ReadonlyState<PD.Values<P>[K]> {
    let field = this.#fields.get(name);
    if (!field) {
      field = derive(this.#values, (v) => v[name]);
      this.#fields.set(name, field);
    }
    return field as ReadonlyState<PD.Values<P>[K]>;
  }

  /** Is the field shown now (its `hideIf` on the current values)? A constant `true` for a field without one. Cached. */
  visible(name: keyof P & string): ReadonlyState<boolean> {
    let visible = this.#visible.get(name);
    if (!visible) {
      const hideIf = this.params[name].hideIf;
      visible = hideIf ? derive(this.#values, (v) => !hideIf(v)) : Always;
      this.#visible.set(name, visible);
    }
    return visible;
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
  fold(path: string, initial = false): State<boolean> {
    let fold = this.#folds.get(path);
    if (!fold) {
      fold = new State(initial);
      this.#folds.set(path, fold);
    }
    return fold;
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
