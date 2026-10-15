/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { Json } from '@molstar/ai-agent';
import type { FieldShape, ObjectShape } from '@molstar/ai-tools';
import { Color } from '@molstar/core/util/color/color';
import { ColorNames } from '@molstar/core/util/color/names';
import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import { camelCaseToWords } from '@molstar/core/util/string';

/**
 * A command's params as a tool's input (agent.md §2; .v6/ai/tools-and-session.md §2): the field shapes the model
 * fills in, and the way back from its arguments to the params' values. Every field may be left out (a param has a
 * default); what JSON cannot carry (files, line graphs, opaque values) is left out of the shape and keeps its default.
 */
export interface ParamsShape {
  readonly shape: ObjectShape;
  /** The params JSON cannot carry: they keep their defaults. */
  readonly omitted: readonly string[];
  /** The arguments (normalized and checked against `shape`) as the params' values, over the defaults. */
  toValues(args: Readonly<Record<string, unknown>>): PD.Values;
  /** The params' values as arguments (what a filled command reads as, for the model). */
  toArgs(values: PD.Values): Record<string, Json>;
  /** Why the arguments cannot be taken (a long select's unknown name, which no enum checks): its choices, all. */
  check(args: Readonly<Record<string, unknown>>): string | undefined;
}

export interface ParamsShapeOptions {
  /** Formats by param name (a command's field hints: "selection", "color"), said in the field's description. */
  readonly formats?: Readonly<Record<string, string | undefined>>;
  /** Options listed by name in a select's description, at most (the enum always has them all). */
  readonly listedOptions?: number;
}

/** Up to this many of a mapped choice's options are named with their labels in the description. */
const ListedOptions = 12;
/**
 * Params every target tool has, said once in the guidance (`ViewerGuidance`) instead of in each tool's schema: what a
 * target is, the selection language, which structure.
 */
const SaidOnce: Readonly<Record<string, (description: string | undefined) => boolean>> = {
  target: (d) => !d || /^What to [a-z ]+: ligand, A:120-130, resn HEM, …$/.test(d),
  structure: (d) => d === 'A structure by label or state ref; empty: the one in view',
  language: () => true,
};
/**
 * A select with more options than this has no enum (a color theme: 47 of them, in every request): its description
 * names a few, and an unknown name fails with all of them (`check`), so the model reads the list only when it needs it.
 */
const LongSelect = 16;
const isLong = (param: PD.Any): param is PD.Select<unknown> =>
  param.type === 'select' && param.options.length > LongSelect;

/** The params as a tool's input. */
export function paramsShape(params: PD.Params, options: ParamsShapeOptions = {}): ParamsShape {
  const fields: Record<string, FieldShape> = {};
  const omitted: string[] = [];
  for (const [name, param] of Object.entries(params)) {
    const field = fieldOf(name, param, options.formats?.[name], options);
    if (field) fields[name] = { ...field, optional: true };
    else omitted.push(name);
  }
  return {
    shape: { kind: 'object', fields },
    omitted,
    toValues: (args) => valuesOf(params, args),
    toArgs: (values) => argsOf(params, values),
    check: (args) => {
      for (const [name, param] of Object.entries(params)) {
        const arg = args[name];
        if (arg === undefined || !isLong(param) || param.options.some(([v]) => String(v) === String(arg))) continue;
        const all = param.options.map(([v]) => String(v)).filter(Boolean);
        return `Unknown ${camelCaseToWords(name).toLowerCase()} '${String(arg)}'. One of: ${all.join(', ')}.`;
      }
      return undefined;
    },
  };
}

// ------------------------------------------------------------------ params → field shapes

/** Whether a param can be given as JSON at all. */
function carried(param: PD.Any): boolean {
  if (param.isHidden) return false;
  switch (param.type) {
    case 'file':
    case 'file-list':
    case 'line-graph':
    case 'value':
      return false;
    case 'converted':
      return carried(param.converted);
    default:
      return true;
  }
}

function fieldOf(
  name: string,
  param: PD.Any,
  format: string | undefined,
  options: ParamsShapeOptions,
): FieldShape | undefined {
  if (!carried(param)) return undefined;
  const shape = shapeOf(param, options);
  if (!shape) return undefined;
  const description = describe(name, param, format, options);
  const fallback = scalarDefault(param);
  return {
    ...shape,
    ...(format && shape.kind === 'string' ? { format } : {}),
    ...(description ? { description } : {}),
    ...(fallback !== undefined ? { default: fallback } : {}),
  };
}

function shapeOf(param: PD.Any, options: ParamsShapeOptions): FieldShape | undefined {
  switch (param.type) {
    case 'select':
      return isLong(param)
        ? { kind: 'string' }
        : { kind: 'string', enum: param.options.map(([value]) => String(value)) };
    case 'multi-select':
      return { kind: 'array', item: { kind: 'string', enum: param.options.map(([value]) => value) } };
    case 'boolean':
      return { kind: 'boolean' };
    case 'text':
      return { kind: 'string', ...(param.multiline ? { verbatim: true } : {}) };
    case 'number':
      return {
        kind: isInteger(param) ? 'integer' : 'number',
        ...(param.min !== undefined ? { min: param.min } : {}),
        ...(param.max !== undefined ? { max: param.max } : {}),
      };
    case 'interval':
      return { kind: 'array', item: { kind: 'number' }, max: 2 };
    case 'vec3':
      return { kind: 'array', item: { kind: 'number' }, max: 3 };
    case 'mat4':
      return { kind: 'array', item: { kind: 'number' }, max: 16 };
    case 'color':
      return { kind: 'string', format: 'color' };
    case 'color-list':
      return {
        kind: 'object',
        fields: {
          kind: { kind: 'string', enum: ['interpolate', 'set'], optional: true },
          colors: { kind: 'array', item: { kind: 'string', format: 'color' } },
        },
      };
    case 'url':
      return { kind: 'string', format: 'url' };
    case 'group':
      return objectOf(param.params, options);
    case 'mapped': {
      const variants: Record<string, ObjectShape> = {};
      for (const [key] of param.select.options) {
        const inner = param.map(key);
        const shape = carried(inner) ? shapeOf(inner, options) : undefined;
        variants[key] = { kind: 'object', fields: shape ? { params: { ...shape, optional: true } } : {} };
      }
      return { kind: 'variant', tag: 'name', variants };
    }
    case 'conditioned': {
      const variants: Record<string, ObjectShape> = {};
      for (const [key, inner] of Object.entries(param.conditionParams as Record<string, PD.Any>)) {
        const shape = carried(inner) ? shapeOf(inner, options) : undefined;
        variants[key] = { kind: 'object', fields: shape ? { value: { ...shape, optional: true } } : {} };
      }
      return { kind: 'variant', tag: 'condition', variants };
    }
    case 'converted':
      return shapeOf(param.converted, options);
    case 'object-list':
      return { kind: 'array', item: objectOf(param.element, options) };
    case 'value-ref':
    case 'data-ref':
      return { kind: 'string', format: 'state node ref' };
    case 'script':
      return {
        kind: 'object',
        fields: {
          language: { kind: 'string', enum: ['pymol', 'mol-script', 'vmd', 'jmol'], optional: true },
          expression: { kind: 'string', verbatim: true },
        },
      };
    default:
      return undefined;
  }
}

function objectOf(params: PD.Params, options: ParamsShapeOptions): ObjectShape {
  const fields: Record<string, FieldShape> = {};
  for (const [name, param] of Object.entries(params)) {
    const field = fieldOf(name, param, undefined, options);
    if (field) fields[name] = { ...field, optional: true };
  }
  return { kind: 'object', fields };
}

const isInteger = (p: PD.Numeric) =>
  p.step !== undefined && Number.isInteger(p.step) && p.step >= 1 && Number.isInteger(p.defaultValue);

/** The default where it is one plain value (scalars and enums); a structure's default is left to its fields. */
function scalarDefault(param: PD.Any): Json | undefined {
  switch (param.type) {
    case 'select':
      return String(param.defaultValue);
    case 'boolean':
    case 'number':
      return param.defaultValue;
    case 'text':
      return param.defaultValue || undefined;
    case 'color':
      return colorText(param.defaultValue);
    default:
      return undefined;
  }
}

/** What the field is, in words: its description or label, and what a select's options mean. */
function describe(
  name: string,
  param: PD.Any,
  format: string | undefined,
  options: ParamsShapeOptions,
): string | undefined {
  const words = camelCaseToWords(name).toLowerCase();
  const label = param.label && param.label.toLowerCase() !== words ? param.label : undefined;
  // A shared param's own words (a structure that may be "every one") stay.
  if (SaidOnce[name]?.(param.description)) return undefined;
  const parts = [param.description ?? label];
  // A select's enum says its names, and a wrong one answers with them: no list in words (every request pays for it).
  if (isLong(param)) parts.push(examplesText(param.options, !!param.description?.includes(',')));
  if (param.type === 'color' || format === 'color') parts.push('A color: a name (red) or hex (#ff0000).');
  if (param.type === 'interval') parts.push('[min, max].');
  if (param.type === 'vec3') parts.push('[x, y, z].');
  if (param.type === 'mat4') parts.push('16 numbers, column-major.');
  if (param.type === 'mapped')
    parts.push(
      `Pick by "name"; its "params" fill in that choice.`,
      optionsText(param.select.options, options.listedOptions ?? ListedOptions),
    );
  if (param.type === 'conditioned') parts.push(`Pick by "condition"; its "value" fills in that choice.`);
  const text = parts.filter(Boolean).join(' ').trim();
  return text || undefined;
}

/** A long select's few first names: an example of the names, not the list (`check` gives it on a wrong one). */
function examplesText(options: PD.Select<unknown>['options'], described: boolean): string {
  const names = options.map(([value]) => String(value)).filter(Boolean);
  const all = `${names.length} names; a wrong one answers with the list`;
  // A description with its own examples (the common ones) says enough of them.
  return described ? `(${all}.)` : `A name: ${names.slice(0, 6).join(', ')}, … (${all}).`;
}

/** A select's options with their labels, where a label says more than the value. */
function optionsText(options: PD.Select<unknown>['options'], max: number): string {
  const named = options.filter(([value, label]) => label && label.toLowerCase() !== String(value).toLowerCase());
  if (!named.length) return '';
  const listed = named.slice(0, max).map(([value, label]) => `${String(value)} (${label})`);
  return `Options: ${listed.join(', ')}${named.length > max ? `, … ${named.length - max} more` : ''}.`;
}

// ------------------------------------------------------------------ arguments → values

const isRecord = (v: unknown): v is Record<string, unknown> => typeof v === 'object' && v !== null && !Array.isArray(v);

function valuesOf(params: PD.Params, args: Readonly<Record<string, unknown>>): PD.Values {
  const values: Record<string, unknown> = {};
  for (const [name, param] of Object.entries(params)) {
    const arg = args[name];
    values[name] = arg === undefined || !carried(param) ? param.defaultValue : valueOf(param, arg);
  }
  return values as PD.Values;
}

function valueOf(param: PD.Any, arg: unknown): unknown {
  switch (param.type) {
    case 'select': {
      const found = param.options.find(([value]) => String(value) === String(arg));
      return found ? found[0] : param.defaultValue;
    }
    case 'color':
      return colorOf(arg) ?? param.defaultValue;
    case 'color-list': {
      if (!isRecord(arg)) return param.defaultValue;
      const colors = Array.isArray(arg.colors) ? arg.colors.map(colorOf).filter((c) => c !== undefined) : [];
      return { kind: arg.kind === 'set' ? 'set' : 'interpolate', colors };
    }
    case 'interval':
    case 'vec3':
    case 'mat4':
      return Array.isArray(arg) ? arg.map(Number) : param.defaultValue;
    case 'group':
      return { ...param.defaultValue, ...valuesOf(param.params, isRecord(arg) ? arg : {}) };
    case 'mapped': {
      const key = isRecord(arg) && typeof arg.name === 'string' ? arg.name : param.defaultValue.name;
      if (!param.select.options.some(([k]) => k === key)) return param.defaultValue;
      const inner = param.map(key);
      const given = isRecord(arg) ? arg.params : undefined;
      const fallback = key === param.defaultValue.name ? param.defaultValue.params : inner.defaultValue;
      return { name: key, params: given === undefined || !carried(inner) ? fallback : valueOf(inner, given) };
    }
    case 'conditioned': {
      const conditions = param.conditionParams as Record<string, PD.Any>;
      const key = isRecord(arg) && typeof arg.condition === 'string' ? arg.condition : undefined;
      const inner = key === undefined ? undefined : conditions[key];
      if (key === undefined || inner === undefined) return param.defaultValue;
      const given = isRecord(arg) ? arg.value : undefined;
      return given === undefined || !carried(inner)
        ? param.conditionedValue(param.defaultValue, key)
        : valueOf(inner, given);
    }
    case 'converted':
      return param.toValue(valueOf(param.converted, arg));
    case 'object-list':
      return Array.isArray(arg)
        ? arg.map((item) => ({ ...param.ctor(), ...valuesOf(param.element, isRecord(item) ? item : {}) }))
        : param.defaultValue;
    case 'value-ref':
    case 'data-ref':
      return PD.Ref(String(arg));
    case 'script':
      return isRecord(arg) && typeof arg.expression === 'string'
        ? {
            language: (['pymol', 'mol-script', 'vmd', 'jmol'] as const).find((l) => l === arg.language) ?? 'pymol',
            expression: arg.expression,
          }
        : param.defaultValue;
    default:
      return arg;
  }
}

/** A color from a name ("red", "light blue"), a hex text ("#ff0000", "ff0000", "0xff0000") or a number. */
export function colorOf(arg: unknown): Color | undefined {
  if (typeof arg === 'number' && Number.isInteger(arg) && arg >= 0 && arg <= 0xffffff) return Color(arg);
  if (typeof arg !== 'string') return undefined;
  const text = arg.trim().toLowerCase();
  const hex = /^(?:#|0x)?([0-9a-f]{6}|[0-9a-f]{3})$/.exec(text)?.[1];
  if (hex) return Color.fromHexStyle(`#${hex.length === 3 ? [...hex].map((c) => c + c).join('') : hex}`);
  const name = text.replace(/[\s_-]+/g, '');
  const found = Object.entries(ColorNames).find(([key]) => key.toLowerCase() === name);
  return found ? (found[1] as Color) : undefined;
}

// ------------------------------------------------------------------ values → arguments

function argsOf(params: PD.Params, values: PD.Values): Record<string, Json> {
  const args: Record<string, Json> = {};
  for (const [name, param] of Object.entries(params)) {
    if (!carried(param)) continue;
    const arg = argOf(param, values[name]);
    if (arg !== undefined) args[name] = arg;
  }
  return args;
}

function argOf(param: PD.Any, value: unknown): Json | undefined {
  if (value === undefined) return undefined;
  switch (param.type) {
    case 'select':
      return String(value);
    case 'color':
      return colorText(value as Color);
    case 'color-list': {
      const list = value as { kind: 'interpolate' | 'set'; colors: unknown[] };
      return {
        kind: list.kind,
        colors: list.colors.map((c) => colorText(Array.isArray(c) ? (c[0] as Color) : (c as Color))),
      };
    }
    case 'interval':
    case 'vec3':
    case 'mat4':
      return Array.from(value as ArrayLike<number>);
    case 'group':
      return argsOf(param.params, value as PD.Values);
    case 'mapped': {
      const { name, params } = value as PD.NamedParams;
      const inner = param.map(name);
      const arg = carried(inner) ? argOf(inner, params) : undefined;
      return arg === undefined ? { name } : { name, params: arg };
    }
    case 'conditioned': {
      const condition = String(param.conditionForValue(value));
      const inner = (param.conditionParams as Record<string, PD.Any>)[condition];
      const arg = inner && carried(inner) ? argOf(inner, value) : undefined;
      return arg === undefined ? { condition } : { condition, value: arg };
    }
    case 'converted':
      return argOf(param.converted, param.fromValue(value));
    case 'object-list':
      return (value as PD.Values[]).map((item) => argsOf(param.element, item));
    case 'value-ref':
    case 'data-ref':
      return (value as PD.Ref).ref;
    case 'script': {
      const script = value as { language: string; expression: string };
      return { language: script.language, expression: script.expression };
    }
    default:
      return value as Json;
  }
}

/** A color as the model reads and writes it: its name when it has one, else hex. */
function colorText(color: Color): string {
  const named = Object.entries(ColorNames).find(([, c]) => c === color);
  return named ? named[0] : Color.toHexStyle(color);
}
