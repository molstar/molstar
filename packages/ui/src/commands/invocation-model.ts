/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// An invocation's values (commands §3.1–3.4, §3.7): how it starts, its state and the views derived from it, setting any
// param, reshaping, moving between fields, the preview. `CommandInvocation` adds typing and running on top.

import { derive, type ReadonlyState, ReactiveModel, State } from '@molstar/core/reactive';
import { sameShallow } from '@molstar/core/util/object';
import type { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import { argsFor, type FieldValue, type FormFieldView, keptValues, readValue } from '../forms/flat-fields.js';
import type { InvocationInit, InvocationState } from './invocation-state.js';
import {
  type Built,
  build,
  type FieldChoices,
  type InvocationOptions,
  isNumeric,
  listChoices,
  missingFields,
  nextEmpty,
  optionsOf,
  type SentenceField,
  sentenceFields,
} from './sentence.js';
import { readTail } from './tail.js';
import type { CommandArgs, CommandDefinition, CommandFormat, CommandValues } from './types.js';

export const errorMessage = (e: unknown) => (e instanceof Error ? e.message : String(e));

export abstract class InvocationModel<Ctx> extends ReactiveModel {
  readonly command: CommandDefinition<Ctx>;
  protected readonly init: InvocationInit<Ctx>;
  #built: Built;
  #values: Record<string, FieldValue> = {};
  #options: Record<string, unknown> = {};
  readonly #state: State<InvocationState>;
  /** Words of the tail no field took: the subclass types them into the first empty field. */
  protected readonly tailRest: readonly string[];
  /** What the preview shows now: the values, the options and the highlighted choice it was drawn for. */
  #previewKey: readonly unknown[] | undefined;
  #previewDeps: readonly unknown[] = [];

  /** The sentence's fields as they read (commands §3.2). */
  readonly sentence: ReadonlyState<readonly SentenceField[]>;
  /** What the active field takes (commands §3.3); none without an active field. */
  readonly choices: ReadonlyState<FieldChoices | undefined>;
  /** The params the sentence does not hold, their values and what differs (commands §3.4). */
  readonly options: ReadonlyState<InvocationOptions>;

  constructor(init: InvocationInit<Ctx>) {
    super();
    const { command, context } = init;
    this.command = command;
    this.init = init;
    // Defaults first, then a preset's values; those that fit, then the tail: typed words win.
    const seed = { ...command.defaults?.(context.value), ...init.preset };
    let built = build(command, context.value, seed);
    const initial = keptValues(built.fields, seed);
    const tail = readTail({
      fields: built.fields,
      tail: init.tail ?? '',
      initial,
      defaulted: new Set(Object.keys(initial)),
      formats: init.formats,
      ctx: context.value,
    });
    Object.assign(initial, tail.values);
    // The tail set a field that reshapes the params ("from PDB-IHM"): the params the values make.
    if (command.reshapes?.some((k) => tail.values[k] !== undefined)) {
      built = build(command, context.value, initial);
      const kept = keptValues(built.fields, initial);
      for (const k of Object.keys(initial)) if (!(k in kept)) delete initial[k];
    }
    this.#built = built;
    this.#values = initial;
    this.tailRest = tail.rest;
    const fields = built.fields;
    const blocked = this.blockedReason();
    // The first empty field, or the last one when all are filled (Enter runs; a slider shows).
    const field = nextEmpty(fields, initial) ?? fields[fields.length - 1]?.id;
    // Compared shallowly: a patch that changes nothing tells nothing.
    this.#state = new State<InvocationState>(
      {
        status: 'editing',
        params: built.params,
        fields,
        formOnly: built.formOnly,
        action: command.action ?? 'Run',
        values: initial,
        optionValues: {},
        ...(field ? { field } : {}),
        text: '',
        pick: -1,
        proposed: [],
        missing: missingFields(fields, initial),
        canAssist: false,
        ...(blocked ? { blocked } : {}),
      },
      sameShallow,
    );

    // Narrow slices of the state: a view recomputes only when what it reads changed (reactive §4: a chain stops at the
    // first unchanged value).
    const slice = <K extends keyof InvocationState>(k: K) => derive(this.#state, (s) => s[k]);
    const sFields = slice('fields');
    const sValues = slice('values');
    const sField = slice('field');
    const sParams = slice('params');
    this.sentence = derive(
      [sParams, sFields, sValues, sField, slice('proposed'), slice('error')],
      (params, fields, values, field, proposed, error) =>
        sentenceFields(params, fields, values, field, proposed, error),
    );
    this.choices = derive([sFields, sValues, sField, slice('text'), context], (fields, values, field, text, ctx) => {
      const f = fields.find((x) => x.id === field);
      return f ? listChoices(f, text, values[f.id], command.fields?.[f.id], this.formatOf(f), ctx) : undefined;
    });
    this.options = derive([sParams, slice('formOnly'), slice('optionValues')], optionsOf);

    // `blocked` follows `when` both ways; the preview follows the context and the values.
    this.listenPermanently(context, () => {
      const blocked = this.blockedReason();
      if (blocked !== this.#state.value.blocked) this.patch({ blocked });
      this.updatePreview();
    });
    this.listenPermanently(this.#state, () => this.updatePreview());
  }

  get state(): ReadonlyState<InvocationState> {
    return this.#state.asReadonly;
  }

  /** The params as they are now. */
  get params(): PD.Params {
    return this.#built.params;
  }

  /** The sentence's values now. */
  values(): Record<string, FieldValue | undefined> {
    return { ...this.#values };
  }

  /** What a run would get now: the sentence's values, the options set, the params' defaults for the rest. */
  args(): CommandArgs {
    return { ...argsFor(this.#built.params, this.#values), ...this.#options };
  }

  // --- Moving between fields (commands §3.2) ------------------------------------------------------

  /** Makes a field of the sentence the one being typed into. */
  focus(field: string): void {
    if (!this.fieldOf(field) || !this.editing()) return;
    this.patch({ field, text: '', pick: -1, canAssist: false, error: undefined });
  }

  /** Moves to the next (or previous) field, wrapping. */
  step(delta: 1 | -1): void {
    if (!this.editing()) return;
    const { fields, field } = this.#state.value;
    const i = fields.findIndex((f) => f.id === field);
    const next = fields[(i + delta + fields.length) % fields.length];
    if (next) this.focus(next.id);
  }

  /** Backspace in an empty field: clears its value, else goes back a field. False at the very start. */
  back(): boolean {
    const s = this.#state.value;
    if (!this.editing() || s.text) return true;
    if (s.field && s.values[s.field] !== undefined) {
      this.clear(s.field);
      return true;
    }
    const i = s.fields.findIndex((f) => f.id === s.field);
    if (i > 0) {
      this.focus((s.fields[i - 1] as FormFieldView).id);
      return true;
    }
    return false;
  }

  // --- Setting values (commands §3.4) ----------------------------------------------------------------

  /**
   * Sets any param: a field of the sentence (checked; the active field moves on to the next empty one unless `stay`),
   * or an option (a form-only or hidden param: what a run, the preview and `check` get over its default). `undefined`
   * clears it. False when not editing, for an unknown name, or for a value the field does not take.
   */
  set(name: string, value: unknown, options: { proposed?: boolean; stay?: boolean } = {}): boolean {
    if (!this.editing()) return false;
    if (value === undefined) return this.clear(name);
    const f = this.fieldOf(name);
    if (!f) {
      if (!(name in this.#built.params)) return false;
      this.#options = { ...this.#options, [name]: value };
      this.patch({ optionValues: this.#options, error: undefined });
      this.#reshape(name);
      return true;
    }
    const read = readValue(f, value as FieldValue);
    if (read === undefined) return this.fail(name, `Not a valid ${f.label}.`);
    if (this.#values[name] !== read) this.#values = { ...this.#values, [name]: read };
    const s = this.#state.value;
    const values = this.#values;
    const active = s.field === name;
    const proposed = options.proposed
      ? s.proposed.includes(name)
        ? s.proposed
        : [...s.proposed, name]
      : s.proposed.filter((x) => x !== name);
    this.patch({
      values,
      field: options.stay || !active ? s.field : (nextEmpty(s.fields, values, name) ?? name),
      ...(active ? { text: '', pick: -1 } : {}),
      canAssist: false,
      error: undefined,
      proposed,
    });
    this.#reshape(name);
    return true;
  }

  /** Empties a field of the sentence, or puts an option back to its default. */
  clear(name: string): boolean {
    if (!this.editing()) return false;
    if (!this.fieldOf(name)) {
      if (!(name in this.#options)) return name in this.#built.params;
      const { [name]: _, ...rest } = this.#options;
      this.#options = rest;
      this.patch({ optionValues: rest });
      this.#reshape(name, true);
      return true;
    }
    if (this.#values[name] !== undefined) {
      const { [name]: _, ...rest } = this.#values;
      this.#values = rest;
    }
    const s = this.#state.value;
    this.patch({
      values: this.#values,
      proposed: s.proposed.filter((x) => x !== name),
      ...(s.field === name ? { pick: -1 } : {}),
    });
    this.#reshape(name, true);
    return true;
  }

  // --- For the subclass ------------------------------------------------------------------------------

  protected get fields(): readonly FormFieldView[] {
    return this.#built.fields;
  }

  protected fieldOf(id: string | undefined): FormFieldView | undefined {
    return id === undefined ? undefined : this.#built.fields.find((f) => f.id === id);
  }

  protected formatOf(f: FormFieldView): CommandFormat<Ctx> | undefined {
    return f.format ? this.init.formats.get(f.format) : undefined;
  }

  protected blockedReason(): string | undefined {
    try {
      const w = this.command.when?.(this.init.context.value);
      return typeof w === 'string' ? w : undefined;
    } catch (e) {
      return errorMessage(e);
    }
  }

  protected editing(): boolean {
    return this.#state.value.status === 'editing';
  }

  protected fail(field: string, text: string, hint?: string): false {
    this.patch({ error: { field, message: text, ...(hint ? { hint } : {}) } });
    return false;
  }

  /** Sets fields of the state; an `undefined` removes one. Nothing changes once closed, except to close. */
  protected patch(partial: Partial<InvocationState>): void {
    const s = this.#state.value;
    if (s.status === 'closed' && partial.status !== 'closed') return;
    const next = { ...s, ...partial } as InvocationState;
    for (const k in partial) {
      const key = k as keyof InvocationState;
      if (partial[key] === undefined) delete next[key];
    }
    if ('values' in partial || 'fields' in partial) next.missing = missingFields(next.fields, next.values);
    this.#state.set(next);
  }

  /** Clears the preview this drew; `always`, whatever is drawn (a run that is done: its effect replaces any preview). */
  protected clearPreview(always = false): void {
    if (always || this.#previewKey !== undefined) this.init.preview?.clear();
    this.#previewKey = undefined;
  }

  /**
   * The preview is what Enter would run (commands §3.7): the values with the highlighted choice, the options, the same
   * context; drawn once the required fields are filled and `when` allows, cleared otherwise. Only while editing: a run
   * keeps it as it was (the run's own edit must not redraw it).
   */
  protected updatePreview(): void {
    const sink = this.init.preview;
    const s = this.#state.value;
    if (!sink || !this.command.preview || s.status !== 'editing') return;
    let values: CommandValues = this.#values;
    const picked = s.pick >= 0 ? this.choices.value?.items[s.pick] : undefined;
    const f = this.fieldOf(s.field);
    if (picked && f && !picked.fill) {
      const v = readValue(f, isNumeric(f) ? Number(picked.value) : picked.value);
      if (v !== undefined) values = { ...values, [f.id]: v };
    }
    const ready = !s.blocked && this.#built.fields.every((x) => !x.required || values[x.id] !== undefined);
    const key = ready ? [this.#values, this.#options, picked?.value] : undefined;
    const ctx = this.init.context.value;
    const deps = this.init.previewDeps?.(ctx) ?? [ctx];
    if (sameItems(key, this.#previewKey) && sameItems(deps, this.#previewDeps)) return;
    const shown = this.#previewKey !== undefined;
    this.#previewKey = key;
    this.#previewDeps = deps;
    if (!key) {
      if (shown) sink.clear();
      return;
    }
    try {
      sink.set(this.command.preview({ ...argsFor(this.#built.params, values), ...this.#options }, ctx));
    } catch {
      // A preview that fails shows nothing.
      this.#previewKey = undefined;
      if (shown) sink.clear();
    }
  }

  /** A param that reshapes the params changed: the params the values make now, keeping what fits. */
  #reshape(changed: string, cleared = false): void {
    if (!this.command.reshapes?.includes(changed)) return;
    let built: Built;
    try {
      built = build(this.command, this.init.context.value, this.#primitiveValues());
    } catch (e) {
      // The params as they were stay: only a field's error says what went wrong.
      this.fail(changed, errorMessage(e));
      return;
    }
    if (built.key === this.#built.key) return;
    this.#built = built;
    this.#values = keptValues(built.fields, this.#values);
    this.#options = Object.fromEntries(Object.entries(this.#options).filter(([k]) => k in built.params));
    // Cleared, it stays on the field (as Backspace does elsewhere); set, it moves on to what is empty.
    const field = cleared ? changed : (nextEmpty(built.fields, this.#values, changed) ?? changed);
    this.patch({
      params: built.params,
      fields: built.fields,
      formOnly: built.formOnly,
      values: this.#values,
      optionValues: this.#options,
      field: built.fields.some((f) => f.id === field) ? field : built.fields[0]?.id,
      text: '',
      pick: -1,
    });
  }

  /** What `params` is built from: the sentence's values and the options that are primitive. */
  #primitiveValues(): CommandValues {
    const out: Record<string, FieldValue> = { ...this.#values };
    for (const [k, v] of Object.entries(this.#options))
      if (typeof v === 'string' || typeof v === 'number' || typeof v === 'boolean') out[k] = v;
    return out;
  }
}

function sameItems(a: readonly unknown[] | undefined, b: readonly unknown[] | undefined): boolean {
  if (a === b) return true;
  if (!a || !b || a.length !== b.length) return false;
  for (let i = 0; i < a.length; i++) if (!Object.is(a[i], b[i])) return false;
  return true;
}
