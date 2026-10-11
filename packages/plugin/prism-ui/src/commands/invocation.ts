/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// One open command: its fields' values, typing, completions, preview, run.
import type { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import {
  argsFor,
  type FieldValue,
  type FormFieldView,
  fits,
  formFields,
  formOnlyParams,
  keptValues,
  readValue,
} from '../forms/flat-fields.js';
import { observableState } from '../reactive/observable-state.js';
import type { ReadonlyState } from '../reactive/readonly-state.js';
import { rankBy } from './match.js';
import { readTail } from './tail.js';
import type {
  CommandArgs,
  CommandAssistant,
  CommandDefinition,
  CommandFormats,
  CommandOutcome,
  CommandResult,
  CommandValues,
  FormatCompletion,
  PreviewSink,
} from './types.js';

export type InvocationStatus = 'editing' | 'resolving' | 'running' | 'cancelling' | 'done' | 'failed' | 'closed';

export interface InvocationState {
  status: InvocationStatus;
  fields: readonly FormFieldView[];
  /** Params edited in the full form only (groups, mapped params, colors, …). */
  formOnly: readonly string[];
  /** The run action's label: "Apply". */
  action: string;
  values: CommandValues;
  /** The field being typed into. */
  field?: string;
  /** What is typed into it (not yet a value). */
  text: string;
  completions: readonly FormatCompletion[];
  /** The highlighted completion; -1 for none. */
  pick: number;
  /** Fields whose values the assistant proposed and nobody has touched since. */
  proposed: readonly string[];
  /** Required fields still empty. */
  missing: readonly string[];
  error?: { field?: string; message: string; hint?: string };
  /** The typed text is not this field's grammar, and an assistant could read it. */
  canAssist: boolean;
  /** What the assistant said about its proposal. */
  note?: string;
  /** The assistant could not express the words as this command. */
  escalated?: string;
  /** The command cannot run now (its `when`). */
  blocked?: string;
  result?: CommandResult;
}

export interface InvocationOptions<Ctx> {
  command: CommandDefinition<Ctx>;
  context: ReadonlyState<Ctx>;
  formats: CommandFormats<Ctx>;
  preview?: PreviewSink;
  assistant?: CommandAssistant<Ctx>;
  /** Text typed after the command's alias: "sticks red". */
  tail?: string;
  /** Values it opens with, over the defaults (a preset's). */
  preset?: CommandValues;
  /** What in the context the preview reads; compared item by item. */
  previewDeps?: (ctx: Ctx) => readonly unknown[];
  /** The run started (its effect may land from now on). */
  onRun?: () => void;
  /** How it ended: once per started run, and for an answer from the assistant. */
  onSettled?: (outcome: CommandOutcome<Ctx>) => void;
}

/** Options first that the text names exactly, by value or label. */
function exactFirst(all: readonly FormatCompletion[], ranked: FormatCompletion[], text: string): FormatCompletion[] {
  const w = text.trim().toLowerCase();
  if (!w) return ranked;
  const exact = all.find((o) => o.value.toLowerCase() === w || o.label.toLowerCase() === w);
  return exact ? [exact, ...ranked.filter((o) => o.value !== exact.value)] : ranked;
}

/** The params as they are now: built again when a field that reshapes them changes. */
interface Built {
  params: PD.Params;
  fields: readonly FormFieldView[];
  /** Identifies the fields' shape, to tell whether a rebuild changed anything. */
  key: string;
}

function build<Ctx>(command: CommandDefinition<Ctx>, ctx: Ctx, values: CommandValues): Built {
  const params = command.params?.(ctx, values);
  if (!params) throw new Error(`Command "${command.id}" has no params.`);
  const fields = formFields(params, command.fields);
  return { params, fields, key: JSON.stringify([fields, formOnlyParams(params, fields)]) };
}

export class CommandInvocation<Ctx> {
  readonly command: CommandDefinition<Ctx>;
  readonly #options: InvocationOptions<Ctx>;
  #built: Built;
  #values: Record<string, FieldValue> = {};
  readonly #state;
  readonly #subs: { unsubscribe(): void }[] = [];
  /** Stops a running assist. */
  #abort = new AbortController();
  /** The run's signal: aborted by an explicit cancel only. */
  #run: AbortController | undefined;
  #previewKey: string | undefined;
  #previewDeps: readonly unknown[] = [];

  constructor(options: InvocationOptions<Ctx>) {
    const { command, context } = options;
    // Defaults first, then a preset's values; those that fit, then the tail: typed words win.
    const seed = { ...command.defaults?.(context.value), ...options.preset };
    let built = build(command, context.value, seed);
    this.command = command;
    this.#options = options;
    const initial = keptValues(built.fields, seed);
    const tail = readTail({
      fields: built.fields,
      tail: options.tail ?? '',
      initial,
      defaulted: new Set(Object.keys(initial)),
      formats: options.formats,
      ctx: context.value,
    });
    Object.assign(initial, tail.values);
    // The tail set a field that reshapes the params ("from PDB-IHM"): the params the values make.
    if (command.reshapes?.some((k) => tail.values[k] !== undefined)) {
      built = build(command, context.value, initial);
      const kept = keptValues(built.fields, initial);
      for (const k of Object.keys(initial)) if (!(k in kept)) delete initial[k];
    }
    const rest = tail.rest;
    this.#built = built;
    this.#values = initial;

    const values = this.values();
    const fields = built.fields;
    const blocked = this.#blockedReason();
    const first = this.#nextEmpty(values);
    const active = first ?? fields[fields.length - 1]?.id;
    this.#state = observableState<InvocationState>({
      status: 'editing',
      fields,
      formOnly: formOnlyParams(built.params, fields),
      action: command.action ?? 'Run',
      values,
      // The first empty field, or the last one when all are filled (Enter runs; a slider shows).
      ...(first ? { field: first } : fields.length ? { field: (fields[fields.length - 1] as FormFieldView).id } : {}),
      text: '',
      completions: this.#listed(active, values),
      pick: -1,
      proposed: [],
      missing: this.#missing(values),
      canAssist: false,
      ...(blocked ? { blocked } : {}),
    });
    // What remains of the tail is typed into the first empty field (and taken when it reads).
    if (rest.length) {
      this.type(rest.join(' '));
      const s = this.#state.asReadonly.value;
      const f = fields.find((x) => x.id === s.field);
      if (f?.format && this.#format(f)?.check(s.text, context.value).ok) this.accept();
    }
    // `blocked` follows `when` both ways; the preview follows the context.
    this.#subs.push(
      context.subscribe(() => {
        const blocked = this.#blockedReason();
        if (blocked !== this.#state.asReadonly.value.blocked) this.#patch({ blocked });
        else this.#updatePreview();
      }),
    );
    this.#updatePreview();
  }

  /** The params as they are now (they follow the values of fields that reshape them). */
  get params(): PD.Params {
    return this.#built.params;
  }

  get state(): ReadonlyState<InvocationState> {
    return this.#state.asReadonly;
  }

  /** The current inline values. */
  values(): Record<string, FieldValue | undefined> {
    return { ...this.#values };
  }

  /** What a run would get now: the inline values plus the params' defaults. */
  args(): CommandArgs {
    return argsFor(this.#built.params, this.#values);
  }

  // --- Editing ---------------------------------------------------------------------------------

  /** Makes a field the one being typed into. */
  focus(field: string): void {
    if (!this.#field(field) || !this.#editing()) return;
    this.#patch({
      field,
      text: '',
      completions: this.#listed(field, this.values()),
      pick: -1,
      canAssist: false,
      error: undefined,
    });
  }

  /** Text typed into the active field: completions follow; nothing is set until accepted. */
  type(text: string): void {
    const s = this.#state.asReadonly.value;
    const f = this.#field(s.field);
    if (!f || !this.#editing()) return;
    const completions = this.#completions(f, text);
    const readable = this.#readable(f, text);
    this.#patch({
      text,
      completions,
      pick: text.trim() && completions.length ? 0 : -1,
      canAssist: !!text.trim() && !completions.length && !readable && this.#assistable(),
      error: undefined,
      escalated: undefined,
      proposed: s.proposed.filter((x) => x !== f.id),
    });
  }

  /** Moves the highlight through the completions; opens them first when a filled field shows none. */
  move(delta: number): void {
    const s = this.#state.asReadonly.value;
    const f = this.#field(s.field);
    if (!s.completions.length && f && s.status === 'editing') {
      const completions = this.#completions(f, s.text);
      if (completions.length) this.#patch({ completions, pick: delta > 0 ? 0 : completions.length - 1 });
      return;
    }
    const n = s.completions.length;
    if (!n) return;
    const pick = s.pick < 0 ? (delta > 0 ? 0 : n - 1) : (((s.pick + delta) % n) + n) % n;
    this.#patch({ pick });
  }

  /** Highlights a completion (e.g. under the pointer). */
  point(index: number): void {
    const s = this.#state.asReadonly.value;
    if (index >= -1 && index < s.completions.length && index !== s.pick) this.#patch({ pick: index });
  }

  /** Takes the highlighted completion, or the typed text when it reads as this field's value. */
  accept(pick?: number): boolean {
    const s = this.#state.asReadonly.value;
    const f = this.#field(s.field);
    if (!f || !this.#editing()) return false;
    const chosen = s.completions[pick ?? s.pick];
    if (chosen?.fill) {
      // An example to edit: into the text, not the value.
      this.type(chosen.value);
      return false;
    }
    if (chosen) return this.set(f.id, isNumeric(f) ? Number(chosen.value) : chosen.value);
    const text = s.text.trim();
    if (!text) return false;
    if (isNumeric(f)) {
      const n = Number(text.replace(/\s*%$/, ''));
      const v = /%$/.test(text) ? n / 100 : n;
      if (!Number.isFinite(v)) return this.#fail(f.id, `"${text}" is not a number.`);
      return this.set(f.id, v);
    }
    if (f.options) {
      const o = f.options.find(
        (x) => x.value.toLowerCase() === text.toLowerCase() || x.label.toLowerCase() === text.toLowerCase(),
      );
      if (o) return this.set(f.id, o.value);
      return this.#fail(
        f.id,
        `"${text}" is not one of the choices.`,
        `Choose: ${f.options.map((x) => x.label).join(', ')}.`,
      );
    }
    const format = this.#format(f);
    if (format) {
      const r = format.check(text, this.#options.context.value);
      if (r.ok) return this.set(f.id, r.value);
      this.#patch({
        error: { field: f.id, message: r.message, ...(r.hint ? { hint: r.hint } : {}) },
        canAssist: !!r.unreadable && this.#assistable(),
      });
      return false;
    }
    return this.set(f.id, text);
  }

  /** Sets a field's value directly (a slider, a pick in 3D); moves on to the next empty field. */
  set(field: string, value: FieldValue, options: { proposed?: boolean; stay?: boolean } = {}): boolean {
    const f = this.#field(field);
    if (!f || !this.#editing()) return false;
    const read = readValue(f, value);
    if (read === undefined) return this.#fail(field, `Not a valid ${f.label}.`);
    this.#values = { ...this.#values, [field]: read };
    const s = this.#state.asReadonly.value;
    const values = this.values();
    const next = options.stay || s.field !== field ? s.field : (this.#nextEmpty(values, field) ?? field);
    this.#patch({
      values,
      field: next,
      text: s.field === field ? '' : s.text,
      completions: s.field === field ? this.#listed(next, values) : s.completions,
      pick: s.field === field ? -1 : s.pick,
      canAssist: false,
      error: undefined,
      proposed: options.proposed ? [...new Set([...s.proposed, field])] : s.proposed.filter((x) => x !== field),
    });
    this.#reshape(field);
    return true;
  }

  /** Empties a field. */
  clear(field: string): void {
    if (!this.#field(field) || !this.#editing()) return;
    const { [field]: _, ...rest } = this.#values;
    this.#values = rest;
    const s = this.#state.asReadonly.value;
    const values = this.values();
    this.#patch({
      values,
      proposed: s.proposed.filter((x) => x !== field),
      ...(s.field === field && !s.text ? { completions: this.#listed(field, values), pick: -1 } : {}),
    });
    this.#reshape(field, true);
  }

  /** Moves to the next (or previous) field. */
  step(delta: 1 | -1): void {
    if (!this.#editing()) return;
    const s = this.#state.asReadonly.value;
    const i = s.fields.findIndex((f) => f.id === s.field);
    const next = s.fields[(i + delta + s.fields.length) % s.fields.length];
    if (next) this.focus(next.id);
  }

  /** Backspace in an empty field: clears its value, else goes back a field. False at the very start. */
  back(): boolean {
    const s = this.#state.asReadonly.value;
    if (!this.#editing() || s.text) return true;
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

  /** Enter: take what is typed, else go to what is missing, else run. */
  async enter(): Promise<CommandResult | undefined> {
    const s = this.#state.asReadonly.value;
    if (s.status !== 'editing') return undefined;
    if (s.text.trim() || s.pick >= 0) {
      // Taking the text (or the highlighted option) runs the command too when it leaves no field empty.
      if (!this.accept()) return undefined;
      const values = this.values();
      if (!this.#built.fields.every((f) => values[f.id] !== undefined)) return undefined;
    }
    return this.confirm();
  }

  /** Tab: let the assistant read unreadable text, else take a completion, else the next field. */
  async tab(delta: 1 | -1 = 1): Promise<void> {
    const s = this.#state.asReadonly.value;
    if (s.canAssist) return this.assist();
    if ((s.text.trim() && s.completions.length) || s.pick >= 0) {
      this.accept();
      return;
    }
    this.step(delta);
  }

  // --- Running ---------------------------------------------------------------------------------

  /**
   * Runs the command with the current values. A run that starts reports one outcome: done (the form is done), failed
   * (it stays open with the error), or cancelled.
   */
  async confirm(): Promise<CommandResult | undefined> {
    const s = this.#state.asReadonly.value;
    if (s.status !== 'editing') return undefined;
    if (s.blocked) {
      this.#patch({ error: { message: s.blocked } });
      return undefined;
    }
    if (s.missing.length) {
      this.#patch({ field: s.missing[0], text: '', completions: this.#listed(s.missing[0], s.values), pick: -1 });
      return undefined;
    }
    this.#patch({ status: 'running', error: undefined });
    const run = new AbortController();
    this.#run = run;
    const ctx = this.#options.context.value;
    const problem = this.#prepare(ctx);
    if (problem) {
      this.#run = undefined;
      this.#patch({ status: 'editing', error: problem });
      return undefined;
    }
    let result: CommandResult | undefined;
    let failure: unknown;
    let failed = false;
    try {
      this.#options.onRun?.();
      result = (await this.command.run(this.args(), ctx, run.signal)) ?? {};
    } catch (e) {
      failed = true;
      failure = e;
    }
    this.#run = undefined;
    // Closed by its owner (its screen went): it reports nothing.
    if (this.#state.asReadonly.value.status === 'closed') return undefined;
    if (!failed) {
      // Done even when cancelled: the effect landed, so the person sees it (and its Undo).
      const done = result as CommandResult;
      this.#patch({ status: 'done', result: done });
      this.#options.preview?.clear();
      this.#previewKey = undefined;
      this.#options.onSettled?.({ kind: 'done', command: this.command, result: done });
      return done;
    }
    if (run.signal.aborted) {
      this.#options.onSettled?.({ kind: 'cancelled', command: this.command });
      return undefined;
    }
    const message = failure instanceof Error ? failure.message : String(failure);
    this.#patch({ status: 'editing', error: { message } });
    this.#options.onSettled?.({ kind: 'failed', command: this.command, message, open: true });
    return undefined;
  }

  /** Cancels the run in flight (an explicit act): its signal aborts; the outcome follows when it settles. */
  cancel(): boolean {
    const s = this.#state.asReadonly.value;
    if (s.status !== 'running' || !this.#run) return false;
    this.#run.abort();
    this.#patch({ status: 'cancelling' });
    return true;
  }

  /** Lets the assistant read the typed words into values (proposed; nothing runs). */
  async assist(): Promise<void> {
    const s = this.#state.asReadonly.value;
    const assistant = this.#options.assistant;
    if (!assistant || !this.command.assist || s.status !== 'editing') return;
    this.#patch({ status: 'resolving', canAssist: false, error: undefined, note: undefined, escalated: undefined });
    const signal = this.#abort.signal;
    try {
      const r = await assistant.propose(
        {
          command: this.command,
          params: this.params,
          values: this.values(),
          text: s.text,
          ...(s.field ? { field: s.field } : {}),
          ctx: this.#options.context.value,
        },
        signal,
      );
      if (signal.aborted || this.#state.asReadonly.value.status === 'closed') return;
      if (r.kind === 'answer') {
        const result: CommandResult = { detail: r.text };
        this.#patch({ status: 'done', result });
        this.#options.onSettled?.({ kind: 'done', command: this.command, result });
        return;
      }
      if (r.kind === 'escalate') {
        this.#patch({ status: 'editing', escalated: s.text, error: { message: r.reason } });
        return;
      }
      this.#patch({ status: 'editing', text: '', completions: [], pick: -1 });
      const rejected: string[] = [];
      for (const [k, v] of Object.entries(r.values)) {
        if (v === undefined) continue;
        const f = this.#field(k);
        const checked =
          f?.format && typeof v === 'string' ? this.#format(f)?.check(v, this.#options.context.value) : undefined;
        if (checked && !checked.ok) {
          rejected.push(`${k}: ${checked.message}`);
          continue;
        }
        if (!this.set(k, checked?.ok ? checked.value : v, { proposed: true })) rejected.push(k);
      }
      const values = this.values();
      this.#patch({
        ...(r.note ? { note: r.note } : {}),
        field: this.#nextEmpty(values) ?? this.#state.asReadonly.value.field,
        ...(rejected.length ? { error: { message: `The proposal did not fit: ${rejected.join('; ')}.` } } : {}),
      });
    } catch (e) {
      if (signal.aborted) return;
      this.#patch({ status: 'editing', error: { message: e instanceof Error ? e.message : String(e) } });
    }
  }

  /**
   * Closes the form at once: drops the preview, stops a running assist. A run in flight is aborted and reports nothing
   * (its owner went); hosts go through the manager, which refuses to close a running form and cancels explicitly.
   */
  close(): void {
    if (this.#state.asReadonly.value.status === 'closed') return;
    this.#abort.abort();
    this.#run?.abort();
    this.#patch({ status: 'closed' });
    if (this.#previewKey !== undefined) this.#options.preview?.clear();
    this.#previewKey = undefined;
    for (const s of this.#subs) s.unsubscribe();
  }

  // --- Internals -------------------------------------------------------------------------------

  /** A field that reshapes the params changed: the params the values make now, keeping what fits. */
  #reshape(changed: string, cleared = false): void {
    if (!this.command.reshapes?.includes(changed)) return;
    let built: Built;
    try {
      built = build(this.command, this.#options.context.value, this.values());
    } catch (e) {
      // The params as they were stay: only a field's error says what went wrong.
      this.#fail(changed, e instanceof Error ? e.message : String(e));
      return;
    }
    if (built.key === this.#built.key) return;
    this.#built = built;
    this.#values = keptValues(built.fields, this.#values);
    const next = this.values();
    // Cleared, it stays on the field (as Backspace does elsewhere); set, it moves on to what is empty.
    const field = cleared ? changed : (this.#nextEmpty(next, changed) ?? changed);
    const f = built.fields.find((x) => x.id === field);
    this.#patch({
      fields: built.fields,
      formOnly: formOnlyParams(built.params, built.fields),
      values: next,
      field,
      text: '',
      completions: f && next[f.id] === undefined ? this.#completions(f, '') : [],
      pick: -1,
    });
  }

  #field(id: string | undefined): FormFieldView | undefined {
    return id === undefined ? undefined : this.#built.fields.find((f) => f.id === id);
  }

  #format(f: FormFieldView) {
    return f.format ? this.#options.formats.get(f.format) : undefined;
  }

  /** What an active field lists before anything is typed: all its options while it is empty. */
  #listed(field: string | undefined, values: CommandValues): FormatCompletion[] {
    const f = this.#field(field);
    return f && values[f.id] === undefined ? this.#completions(f, '') : [];
  }

  #completions(f: FormFieldView, text: string): FormatCompletion[] {
    // Options match by what they read as (their values may be opaque ids).
    if (f.options) {
      const all = f.options.map((o) => ({ value: o.value, label: o.label }));
      return exactFirst(
        all,
        rankBy(all, text, (o) => o.label),
        text,
      );
    }
    const format = this.#format(f);
    if (!format) return [];
    const ctx = this.#options.context.value;
    const listed = format.complete(text, ctx);
    const completions = exactFirst(listed, listed, text);
    // Text that reads as a value comes first, as read ("A 120-133, 155" → "A:120-133, A:155"), so Enter takes what
    // was typed rather than a name that happens to match.
    const typed = text.trim() ? format.check(text, ctx) : undefined;
    if (typed?.ok && !completions.some((c) => c.value === typed.value))
      return [{ value: typed.value, label: typed.value, detail: 'as typed' }, ...completions];
    return completions;
  }

  /** Would the text read as a value of this field (then there is nothing for the assistant)? */
  #readable(f: FormFieldView, text: string): boolean {
    if (!text.trim()) return true;
    if (isNumeric(f)) return Number.isFinite(Number(text.replace(/\s*%$/, '')));
    if (f.options) return false;
    const format = this.#format(f);
    if (!format) return true;
    const r = format.check(text, this.#options.context.value);
    return r.ok || !r.unreadable;
  }

  #assistable(): boolean {
    return !!this.command.assist && !!this.#options.assistant;
  }

  #missing(values: CommandValues): string[] {
    return this.#built.fields.filter((f) => f.required && values[f.id] === undefined).map((f) => f.id);
  }

  /** The first empty required field (after `after`, wrapping), then the first empty optional one. */
  #nextEmpty(values: CommandValues, after?: string): string | undefined {
    const fields = this.#built.fields;
    const i = after ? fields.findIndex((f) => f.id === after) : -1;
    const order = [...fields.slice(i + 1), ...fields.slice(0, i + 1)];
    return (
      order.find((f) => f.required && values[f.id] === undefined)?.id ??
      order.find((f) => values[f.id] === undefined && f.id !== after)?.id
    );
  }

  #blockedReason(): string | undefined {
    try {
      const w = this.command.when?.(this.#options.context.value);
      return typeof w === 'string' ? w : undefined;
    } catch (e) {
      return e instanceof Error ? e.message : String(e);
    }
  }

  #editing(): boolean {
    return this.#state.asReadonly.value.status === 'editing';
  }

  #fail(field: string, message: string, hint?: string): false {
    this.#patch({ error: { field, message, ...(hint ? { hint } : {}) } });
    return false;
  }

  /** Why the run cannot start: `when`, a value its format no longer reads, or the command's own check. */
  #prepare(ctx: Ctx): InvocationState['error'] | undefined {
    const w = this.command.when?.(ctx);
    if (typeof w === 'string') return { message: w };
    for (const f of this.#built.fields) {
      const v = this.#values[f.id];
      const format = this.#format(f);
      if (format && typeof v === 'string') {
        const r = format.check(v, ctx);
        if (!r.ok) return { message: r.message, ...(r.hint ? { hint: r.hint } : {}) };
      }
      if (v !== undefined && !fits(f, v)) return { field: f.id, message: `Not a valid ${f.label}.` };
    }
    const problem = this.command.check?.(this.args(), ctx);
    return problem ? { message: problem } : undefined;
  }

  #patch(partial: Partial<InvocationState>): void {
    const s = this.#state.asReadonly.value;
    if (s.status === 'closed' && partial.status !== 'closed') return;
    const next = { ...s, ...partial } as InvocationState;
    for (const k of Object.keys(partial) as (keyof InvocationState)[]) if (partial[k] === undefined) delete next[k];
    if ('values' in partial) next.missing = this.#missing(next.values);
    this.#state.set(next);
    this.#updatePreview();
  }

  #updatePreview(): void {
    const sink = this.#options.preview;
    const s = this.#state?.asReadonly.value;
    if (!sink || !this.command.preview || !s) return;
    // While it runs the preview stays as it was (the run's own edit must not redraw it); done clears it.
    if (s.status === 'closed' || s.status === 'done' || s.status === 'running' || s.status === 'cancelling') return;
    const values: Record<string, FieldValue | undefined> = { ...s.values };
    const picked = s.completions[s.pick];
    const f = this.#field(s.field);
    if (picked && f) values[f.id] = readValue(f, isNumeric(f) ? Number(picked.value) : picked.value);
    const ready = this.#built.fields.every((x) => !x.required || values[x.id] !== undefined);
    const key = ready && !s.blocked ? JSON.stringify(values) : undefined;
    // The preview is what Enter would run: the same values and the same context.
    const ctx = this.#options.context.value;
    const deps = this.#options.previewDeps?.(ctx) ?? [ctx];
    const same =
      key === this.#previewKey &&
      deps.length === this.#previewDeps.length &&
      deps.every((d, i) => Object.is(d, this.#previewDeps[i]));
    if (same) return;
    const shown = this.#previewKey !== undefined;
    this.#previewKey = key;
    this.#previewDeps = deps;
    if (key === undefined) {
      if (shown) sink.clear();
    } else sink.set(this.command.preview(argsFor(this.#built.params, values), ctx));
  }
}

function isNumeric(f: FormFieldView): boolean {
  return f.kind === 'number' || f.kind === 'slider';
}
