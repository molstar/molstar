/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { Json } from '@molstar/ai-agent';
import type { z } from 'zod';
import { type ToolInput, zodBacked } from './input.js';
import type { FieldShape, ObjectShape } from './shape.js';

/**
 * Tools declared with zod (tools-and-session.md §1): `input: zodInput(z.object({ … }))`. The object is the validator
 * and the JSON Schema's source as it is; for normalizing the arguments it is also read as a field shape (§2).
 */

interface Def {
  type: string;
  innerType?: z.ZodType;
  element?: z.ZodType;
  entries?: Record<string, string>;
  values?: unknown[];
  shape?: Record<string, z.ZodType>;
  options?: z.ZodType[];
  discriminator?: string;
  defaultValue?: Json;
}
const def = (s: z.ZodType): Def => (s as unknown as { _zod: { def: Def } })._zod.def;

const VERBATIM = new WeakSet<z.ZodType>();

/** A field taken as sent (code, a snippet to find): not trimmed, and an empty string is a value. */
export function verbatim<T extends z.ZodType>(s: T): T {
  VERBATIM.add(s);
  return s;
}

const WRAPPERS = ['optional', 'nullable', 'default', 'prefault', 'readonly', 'catch'];

const bound = (v: number | null | undefined, safe: number) =>
  v !== null && v !== undefined && Number.isFinite(v) && v !== safe;

function shapeOf(schema: z.ZodType): FieldShape {
  let optional = false;
  let dflt: Json | undefined;
  let description: string | undefined;
  let marked = false;
  let cur = schema;
  for (;;) {
    description ??= cur.description;
    marked ||= VERBATIM.has(cur);
    const d = def(cur);
    if (!WRAPPERS.includes(d.type) || !d.innerType) break;
    if (d.type === 'optional') optional = true;
    if (d.type === 'default' || d.type === 'prefault') dflt = d.defaultValue;
    cur = d.innerType;
  }
  const d = def(cur);
  const extra = {
    ...(optional ? { optional } : {}),
    ...(dflt !== undefined ? { default: dflt } : {}),
    ...(description ? { description } : {}),
  };
  const number = cur as unknown as { isInt?: boolean; minValue?: number | null; maxValue?: number | null };
  switch (d.type) {
    case 'string':
      return { kind: 'string', ...(marked ? { verbatim: true } : {}), ...extra };
    case 'number':
      return {
        kind: number.isInt ? 'integer' : 'number',
        // An integer's safe range is not a bound anyone asked for.
        ...(bound(number.minValue, Number.MIN_SAFE_INTEGER) ? { min: number.minValue as number } : {}),
        ...(bound(number.maxValue, Number.MAX_SAFE_INTEGER) ? { max: number.maxValue as number } : {}),
        ...extra,
      };
    case 'boolean':
      return { kind: 'boolean', ...extra };
    case 'enum':
      return { kind: 'string', enum: Object.values(d.entries ?? {}), ...extra };
    case 'literal':
      if (d.values?.every((v) => typeof v === 'string'))
        return { kind: 'string', enum: d.values as string[], ...extra };
      return d.values?.every((v) => typeof v === 'number') ? { kind: 'number', ...extra } : { kind: 'json', ...extra };
    case 'array':
      return { kind: 'array', item: d.element ? shapeOf(d.element) : { kind: 'json' }, ...extra };
    case 'object':
      return { ...objectShapeOf(d.shape ?? {}), ...extra };
    case 'union': {
      const tag = d.discriminator;
      if (tag === undefined) return { kind: 'json', ...extra };
      const variants: Record<string, ObjectShape> = {};
      for (const option of d.options ?? []) {
        const fields = def(option).shape ?? {};
        const literal = fields[tag] ? def(fields[tag]).values : undefined;
        const { [tag]: _tag, ...rest } = fields;
        for (const key of literal ?? []) variants[String(key)] = objectShapeOf(rest);
      }
      return { kind: 'variant', tag, variants, ...extra };
    }
    default:
      return { kind: 'json', ...extra };
  }
}

function objectShapeOf(fields: Record<string, z.ZodType>): ObjectShape {
  return { kind: 'object', fields: Object.fromEntries(Object.entries(fields).map(([k, f]) => [k, shapeOf(f)])) };
}

/**
 * A tool's input from a zod object: its JSON Schema (draft 2020-12, the input side) and its check are zod's own, its
 * field shape is read from the object so arguments are normalized as for any tool. Mark code-like strings with
 * `verbatim(…)`.
 */
export function zodInput<S extends z.ZodObject>(schema: S): ToolInput<z.output<S>> {
  return zodBacked<z.output<S>>(objectShapeOf(schema.shape as Record<string, z.ZodType>), schema);
}
