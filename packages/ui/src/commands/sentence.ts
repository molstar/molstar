/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// What an invocation shows, computed from its values (commands §3.2–3.4): its params split into the sentence and the
// options, each field as it reads, the active field's choices. Pure functions; the invocation derives its states from
// them.

import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import { camelCaseToWords } from '@molstar/core/util/string';
import {
  type FieldKind,
  type FieldValue,
  type FormFieldView,
  fits,
  formFields,
  formOnlyParams,
  readValue,
  valueLabel,
} from '../forms/flat-fields.js';
import { rankBy } from './match.js';
import { checkFormat } from './format.js';
import type { CommandDefinition, CommandFieldHint, CommandFormat, CommandValues, FormatCompletion } from './types.js';

/** The params as they are now: built again when a field that reshapes them changes. */
export interface Built {
  params: PD.Params;
  /** The sentence's fields, in order. */
  fields: readonly FormFieldView[];
  /** The params the sentence does not hold (hidden ones not), in order. */
  formOnly: readonly string[];
  /** Identifies the shape, to tell whether a rebuild changed anything. */
  key: string;
}

export function build<Ctx>(command: CommandDefinition<Ctx>, ctx: Ctx, values: CommandValues): Built {
  const params = command.params?.(ctx, values);
  if (!params) throw new Error(`Command "${command.id}" has no params.`);
  const fields = formFields(params, command.fields);
  const formOnly = formOnlyParams(params, fields);
  return { params, fields, formOnly, key: JSON.stringify([fields, formOnly]) };
}

export function isNumeric(f: FormFieldView): boolean {
  return f.kind === 'number' || f.kind === 'slider';
}

/** Typed text as a number: "0.4", and "60%" as 0.6. */
export function readNumber(text: string): number | undefined {
  const t = text.trim();
  const n = Number(t.replace(/\s*%$/, ''));
  if (!t || !Number.isFinite(n)) return undefined;
  return /%$/.test(t) ? n / 100 : n;
}

export function paramLabel(params: PD.Params, name: string): string {
  return params[name]?.label || camelCaseToWords(name);
}

/** Required fields still empty, in order. */
export function missingFields(fields: readonly FormFieldView[], values: CommandValues): string[] {
  const out: string[] = [];
  for (const f of fields) if (f.required && values[f.id] === undefined) out.push(f.id);
  return out;
}

/** The first empty required field (after `after`, wrapping), then the first empty optional one. */
export function nextEmpty(fields: readonly FormFieldView[], values: CommandValues, after?: string): string | undefined {
  const i = after ? fields.findIndex((f) => f.id === after) : -1;
  const n = fields.length;
  let optional: string | undefined;
  for (let k = 1; k <= n; k++) {
    const f = fields[(i + k) % n] as FormFieldView;
    if (values[f.id] !== undefined || f.id === after) continue;
    if (f.required) return f.id;
    optional ??= f.id;
  }
  return optional;
}

// --- The sentence (commands §3.2) ---------------------------------------------------------------

/** A field of the sentence as the view draws it. */
export interface SentenceField extends FormFieldView {
  /** Its value, if set. */
  value: FieldValue | undefined;
  /** The value as it reads (an option's label, on/off, the number); '' when empty. */
  display: string;
  /** Its param's default (what an empty optional field stands for); none for a required text. */
  defaultValue: FieldValue | undefined;
  defaultDisplay: string;
  empty: boolean;
  /** The field being typed into. */
  active: boolean;
  /** The assistant proposed its value and nobody touched it since. */
  proposed: boolean;
  /** Set to something other than its param's default. */
  differs: boolean;
  /** The error is this field's. */
  error?: string;
}

export function sentenceFields(
  params: PD.Params,
  fields: readonly FormFieldView[],
  values: CommandValues,
  active: string | undefined,
  proposed: readonly string[],
  error: { field?: string; message: string } | undefined,
): SentenceField[] {
  return fields.map((f) => {
    const value = values[f.id];
    const param = params[f.id];
    const defaultValue = param ? readValue(f, param.defaultValue as FieldValue) : undefined;
    const differs =
      value !== undefined && (defaultValue === undefined || !PD.isParamEqual(param as PD.Any, value, defaultValue));
    return {
      ...f,
      value,
      display: valueLabel(f, value),
      defaultValue,
      defaultDisplay: valueLabel(f, defaultValue),
      empty: value === undefined,
      active: f.id === active,
      proposed: proposed.includes(f.id),
      differs,
      ...(error?.field === f.id ? { error: error.message } : {}),
    };
  });
}

// --- The options (commands §3.4) ----------------------------------------------------------------

export interface ChangedOption {
  name: string;
  label: string;
  value: unknown;
}

/** The params the sentence does not hold, for a params form (the menu's Options fold). */
export interface InvocationOptions {
  /** The form-only params, in order. */
  params: PD.Params;
  /** Their values now: those set over the defaults. */
  values: Record<string, unknown>;
  /** Those that differ from their default, in order: the fold's summary. */
  changed: readonly ChangedOption[];
}

export function optionsOf(
  params: PD.Params,
  formOnly: readonly string[],
  set: Readonly<Record<string, unknown>>,
): InvocationOptions {
  const own: PD.Params = {};
  for (const name of formOnly) own[name] = params[name];
  const defaults = PD.getDefaultValues(own) as Record<string, unknown>;
  const values: Record<string, unknown> = { ...defaults };
  const changed: ChangedOption[] = [];
  for (const name of formOnly) {
    if (!(name in set)) continue;
    const value = set[name];
    values[name] = value;
    if (!PD.isParamEqual(own[name], value, defaults[name]))
      changed.push({ name, label: paramLabel(params, name), value });
  }
  return { params: own, values, changed };
}

// --- The active field's choices (commands §3.3) ------------------------------------------------------

export interface FieldChoice extends FormatCompletion {
  /** It is the field's value now. */
  current: boolean;
}

/** What the active field takes, for the view to list: always the active field's, filtered by what is typed. */
export interface FieldChoices {
  field: string;
  kind: FieldKind;
  /** What is typed into the field. */
  text: string;
  /**
   * A select's options (an exact name first, then by match), a boolean's two values, a format's completions (what the
   * text reads as first), a number's presets.
   */
  items: readonly FieldChoice[];
  /** A format's hint: what can be typed. */
  hint?: string;
  /** A number's bounds, for a slider. */
  range?: { min?: number; max?: number; step?: number };
}

/** A slider field's bounds: what a view draws a range from (`InvocationModel.range`). */
export interface FieldRange {
  field: string;
  min: number;
  max: number;
  step?: number;
}

/** The field's range when it is a slider with both bounds; none otherwise. */
export function fieldRange(f: FormFieldView | undefined): FieldRange | undefined {
  if (!f || f.kind !== 'slider' || f.min === undefined || f.max === undefined) return undefined;
  return { field: f.id, min: f.min, max: f.max, ...(f.step !== undefined ? { step: f.step } : {}) };
}

/** Options first that the text names exactly, by value or label. */
function exactFirst(all: readonly FormatCompletion[], ranked: FormatCompletion[], text: string): FormatCompletion[] {
  const w = text.trim().toLowerCase();
  if (!w) return ranked;
  const exact = all.find((o) => o.value.toLowerCase() === w || o.label.toLowerCase() === w);
  return exact ? [exact, ...ranked.filter((o) => o.value !== exact.value)] : ranked;
}

function formatCompletions<Ctx>(format: CommandFormat<Ctx>, text: string, ctx: Ctx): FormatCompletion[] {
  let listed: FormatCompletion[];
  try {
    listed = format.complete(text, ctx);
  } catch {
    listed = [];
  }
  const completions = exactFirst(listed, listed, text);
  // Text that reads as a value comes first, as read ("A 120-133, 155" → "A:120-133, A:155"), so Enter takes what was
  // typed rather than a name that happens to match.
  const typed = text.trim() ? checkFormat(format, text, ctx) : undefined;
  if (typed?.ok && !completions.some((c) => c.value === typed.value))
    return [{ value: typed.value, label: typed.value, detail: 'as typed' }, ...completions];
  return completions;
}

export function listChoices<Ctx>(
  f: FormFieldView,
  text: string,
  value: FieldValue | undefined,
  hint: CommandFieldHint | undefined,
  format: CommandFormat<Ctx> | undefined,
  ctx: Ctx,
): FieldChoices {
  const base: FieldChoices = { field: f.id, kind: f.kind, text, items: [] };
  const mark = (items: readonly FormatCompletion[], same: (c: FormatCompletion) => boolean) =>
    items.map((c) => ({ ...c, current: value !== undefined && same(c) }));
  // Options match by what they read as (their values may be opaque ids).
  if (f.options) {
    const all = f.options.map((o) => ({ value: o.value, label: o.label }));
    const items = exactFirst(
      all,
      rankBy(all, text, (o) => o.label),
      text,
    );
    return { ...base, items: mark(items, (c) => c.value === String(value)) };
  }
  if (isNumeric(f)) {
    const presets: FormatCompletion[] = [];
    for (const p of hint?.presets ?? []) {
      const [v, label] = typeof p === 'number' ? [p, String(p)] : p;
      if (fits(f, v)) presets.push({ value: String(v), label });
    }
    const range = {
      ...(f.min !== undefined ? { min: f.min } : {}),
      ...(f.max !== undefined ? { max: f.max } : {}),
      ...(f.step !== undefined ? { step: f.step } : {}),
    };
    return {
      ...base,
      items: mark(
        rankBy(presets, text, (p) => p.label),
        (c) => Number(c.value) === value,
      ),
      ...(f.kind === 'slider' ? { range } : {}),
    };
  }
  if (!format) return base;
  let formatHint: string | undefined;
  try {
    formatHint = format.hint?.(ctx);
  } catch {
    formatHint = undefined;
  }
  return {
    ...base,
    items: mark(formatCompletions(format, text, ctx), (c) => c.value === value),
    ...(formatHint ? { hint: formatHint } : {}),
  };
}

/** Would the text read as a value of this field (then there is nothing for the assistant)? */
export function isReadable<Ctx>(
  f: FormFieldView,
  text: string,
  format: CommandFormat<Ctx> | undefined,
  ctx: Ctx,
): boolean {
  if (!text.trim()) return true;
  if (isNumeric(f)) return readNumber(text) !== undefined;
  if (f.options) return false;
  if (!format) return true;
  const r = checkFormat(format, text, ctx);
  return r.ok || !r.unreadable;
}
