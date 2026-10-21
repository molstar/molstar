/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import { camelCaseToWords } from '@molstar/core/util/string';

/** A value an inline field holds. */
export type FieldValue = string | number | boolean;
export type FieldValues = Readonly<Record<string, FieldValue | undefined>>;

export interface FormOption {
  value: string;
  label: string;
}

/** How a param reads inline, as part of a sentence: "show <what> as <representation>". */
export interface FieldHint {
  /** A quiet word before the field: "as", "in", "onto". */
  before?: string;
  /** A named input type ("selection", "color") that the host checks and completes. */
  format?: string;
}

export type FieldKind = 'text' | 'select' | 'boolean' | 'number' | 'slider';

/** A param that can be edited inline, as hosts need it. */
export interface FormFieldView {
  /** The param name. */
  id: string;
  kind: FieldKind;
  label: string;
  before?: string;
  format?: string;
  /** Closed choices: select options, or on/off for a boolean. */
  options?: readonly FormOption[];
  required: boolean;
  min?: number;
  max?: number;
  step?: number;
}

const BooleanOptions: readonly FormOption[] = [
  { value: 'true', label: 'on' },
  { value: 'false', label: 'off' },
];

/**
 * The params that read as inline fields, in order. Text, string selects, booleans and numbers are inline; everything
 * else (groups, mapped params, colors, …) is edited in the full form only, and hidden params never show.
 */
export function formFields(params: PD.Params, hints?: Readonly<Record<string, FieldHint>>): FormFieldView[] {
  const fields: FormFieldView[] = [];
  for (const id of Object.keys(params)) {
    const param = params[id];
    if (param.isHidden) continue;
    const field = inlineField(id, param);
    if (!field) continue;
    const hint = hints?.[id];
    if (hint?.before) field.before = hint.before;
    if (hint?.format) field.format = hint.format;
    fields.push(field);
  }
  return fields;
}

function inlineField(id: string, param: PD.Any): FormFieldView | undefined {
  const base = { id, label: param.label || camelCaseToWords(id), required: !param.isOptional };
  switch (param.type) {
    case 'text':
      return { ...base, kind: 'text' };
    case 'boolean':
      return { ...base, kind: 'boolean', options: BooleanOptions };
    case 'select':
      if (!param.options.every((o) => typeof o[0] === 'string')) return undefined;
      return {
        ...base,
        kind: 'select',
        options: param.options.map(([value, label]) => ({ value: value as string, label })),
      };
    case 'number': {
      const slider = param.min !== undefined && param.max !== undefined;
      return {
        ...base,
        kind: slider ? 'slider' : 'number',
        ...(param.min !== undefined ? { min: param.min } : {}),
        ...(param.max !== undefined ? { max: param.max } : {}),
        ...(param.step !== undefined ? { step: param.step } : {}),
      };
    }
    default:
      return undefined;
  }
}

/** The names of params that are edited in the full form only (hidden ones are not offered anywhere). */
export function formOnlyParams(params: PD.Params, fields: readonly FormFieldView[]): string[] {
  const inline = new Set(fields.map((f) => f.id));
  return Object.keys(params).filter((k) => !inline.has(k) && !params[k].isHidden);
}

/** Is the field's value a string (so words can name it)? */
export function isStringField(field: FormFieldView): boolean {
  return field.kind === 'text' || field.kind === 'select' || field.kind === 'boolean';
}

/**
 * The value as the field holds it, or undefined when it does not fit: an option of a select, on/off (or a boolean)
 * for a boolean, a finite number within bounds, non-empty text for a required text field.
 */
export function readValue(field: FormFieldView, value: FieldValue | undefined): FieldValue | undefined {
  if (value === undefined) return undefined;
  switch (field.kind) {
    case 'text':
      return typeof value === 'string' && (!field.required || value.length > 0) ? value : undefined;
    case 'select':
      return typeof value === 'string' && field.options?.some((o) => o.value === value) ? value : undefined;
    case 'boolean':
      if (typeof value === 'boolean') return value;
      return value === 'true' ? true : value === 'false' ? false : undefined;
    case 'number':
    case 'slider': {
      if (typeof value !== 'number' || !Number.isFinite(value)) return undefined;
      if (field.min !== undefined && value < field.min) return undefined;
      if (field.max !== undefined && value > field.max) return undefined;
      return value;
    }
  }
}

/** Does the value fit the field? */
export function fits(field: FormFieldView, value: FieldValue | undefined): boolean {
  return readValue(field, value) !== undefined;
}

/** The values the fields take (those that fit), as the fields hold them. */
export function keptValues(fields: readonly FormFieldView[], values: FieldValues): Record<string, FieldValue> {
  const out: Record<string, FieldValue> = {};
  for (const [k, v] of Object.entries(values)) {
    const f = fields.find((x) => x.id === k);
    const read = f ? readValue(f, v) : undefined;
    if (read !== undefined) out[k] = read;
  }
  return out;
}

/** The text a value shows as: an option's label, on/off, or the value itself. */
export function valueLabel(field: FormFieldView, value: FieldValue | undefined): string {
  if (value === undefined) return '';
  const option = field.options?.find((o) => o.value === String(value));
  return option ? option.label : String(value);
}

/**
 * The arguments a run gets: the inline values, then the params' defaults for what is not set inline (form-only
 * params and optional fields left empty stay as PD defines them: optional ones are absent).
 */
export function argsFor(params: PD.Params, values: FieldValues): Record<string, unknown> {
  const args: Record<string, unknown> = { ...PD.getDefaultValues(params) };
  for (const [k, v] of Object.entries(values)) {
    if (v !== undefined) args[k] = v;
  }
  for (const k of Object.keys(params)) {
    if (values[k] === undefined && params[k].isOptional) delete args[k];
  }
  return args;
}
