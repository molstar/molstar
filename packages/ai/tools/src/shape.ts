/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { Json } from '@molstar/ai-agent';

/**
 * A tool's input as a small neutral shape (tools-and-session.md §2): what a host declares without importing a
 * validator, what the arguments are normalized against, and what `ai-tools` makes its JSON Schema and check from.
 */
export type FieldShape = (
  | {
      readonly kind: 'string';
      readonly enum?: readonly string[];
      /** Taken as sent (code, a snippet to find): not trimmed, and an empty string is a value. */
      readonly verbatim?: boolean;
      /** What the text looks like ("A:315, B:200-210"): said in the description, not checked. */
      readonly format?: string;
    }
  | { readonly kind: 'number' | 'integer'; readonly min?: number; readonly max?: number }
  | { readonly kind: 'boolean' }
  | { readonly kind: 'array'; readonly item: FieldShape; readonly max?: number }
  | ObjectShape
  | { readonly kind: 'variant'; readonly tag: string; readonly variants: Readonly<Record<string, ObjectShape>> }
  | { readonly kind: 'json' }
) & { readonly optional?: boolean; readonly description?: string; readonly default?: Json };

export interface ObjectShape {
  readonly kind: 'object';
  readonly fields: Readonly<Record<string, FieldShape>>;
}

// ------------------------------------------------------------------ the arguments a shape gives

type Fields<S extends { readonly fields: Readonly<Record<string, FieldShape>> }> = S['fields'];
type OptionalKeys<F extends Readonly<Record<string, FieldShape>>> = {
  [K in keyof F]: F[K] extends { readonly optional: true } ? K : never;
}[keyof F];

/** The value a field shape checks to (a field with a default is always there). */
export type ValueOf<F extends FieldShape> = F extends { readonly kind: 'string'; readonly enum: readonly (infer E)[] }
  ? E
  : F extends { readonly kind: 'string' }
    ? string
    : F extends { readonly kind: 'number' | 'integer' }
      ? number
      : F extends { readonly kind: 'boolean' }
        ? boolean
        : F extends { readonly kind: 'array'; readonly item: infer I extends FieldShape }
          ? ValueOf<I>[]
          : F extends ObjectShape
            ? ArgsOf<F>
            : F extends {
                  readonly kind: 'variant';
                  readonly tag: infer T extends string;
                  readonly variants: infer V extends Readonly<Record<string, ObjectShape>>;
                }
              ? { [K in keyof V]: { [P in T]: K } & ArgsOf<V[K]> }[keyof V]
              : Json;

/** The arguments an object shape checks to, as `run` gets them. */
export type ArgsOf<S extends ObjectShape> = {
  [K in Exclude<keyof Fields<S>, OptionalKeys<Fields<S>>>]: ValueOf<Fields<S>[K]>;
} & { [K in OptionalKeys<Fields<S>>]?: ValueOf<Fields<S>[K]> };

// ------------------------------------------------------------------ normalizing

const isRecord = (v: unknown): v is Record<string, unknown> => typeof v === 'object' && v !== null && !Array.isArray(v);

const trimmed = (v: unknown) => (typeof v === 'string' ? v.trim() : v);

const sameCase = (values: readonly string[], v: string) => values.find((x) => x.toLowerCase() === v.toLowerCase()) ?? v;

/**
 * What small models send, made what the shape asks for, before it is checked (tools-and-session.md §2): an empty
 * string or null is absent; strings are trimmed (not verbatim ones); "3" is a number and "true" a boolean where one is
 * asked for, a number a string where a string is; a list of strings is one comma-joined text where one text is asked
 * for, and "A, B" a list where a list of non-objects is; an enum's value in another case is its value; objects (in a
 * list too) are read by their own fields, a variant's by the variant its tag names.
 */
export function normalize(shape: ObjectShape, args: Readonly<Record<string, unknown>>): Record<string, unknown> {
  const out: Record<string, unknown> = {};
  for (const [key, raw] of Object.entries(args)) {
    const field = shape.fields[key];
    if (field?.kind === 'string' && field.verbatim && typeof raw === 'string') {
      out[key] = raw;
      continue;
    }
    const v = trimmed(raw);
    if (v === null || v === '') continue;
    out[key] = field ? coerce(field, v) : v;
  }
  return out;
}

function coerce(field: FieldShape, value: unknown): unknown {
  let v = value;
  switch (field.kind) {
    case 'number':
    case 'integer':
      return typeof v === 'string' && Number.isFinite(Number(v)) ? Number(v) : v;
    case 'boolean':
      return v === 'true' || v === 'false' ? v === 'true' : v;
    case 'string':
      if (typeof v === 'number') v = String(v);
      // A list where one text is asked ("residues": ["A:315", "A:318"]) is that text.
      else if (Array.isArray(v) && v.every((x) => typeof x === 'string')) v = v.map((x: string) => x.trim()).join(', ');
      return field.enum && typeof v === 'string' ? sameCase(field.enum, v) : v;
    case 'array': {
      const item = field.item;
      let list = v;
      if (typeof v === 'string' && ['string', 'number', 'integer', 'boolean'].includes(item.kind))
        list = v
          .split(',')
          .map((x) => x.trim())
          .filter(Boolean);
      return Array.isArray(list) ? list.map((x) => coerceItem(item, x)) : list;
    }
    case 'object':
      return isRecord(v) ? normalize(field, v) : v;
    case 'variant': {
      if (!isRecord(v)) return v;
      const tag = trimmed(v[field.tag]);
      const key =
        typeof tag === 'string'
          ? Object.keys(field.variants).find((k) => k.toLowerCase() === tag.toLowerCase())
          : undefined;
      const variant = key === undefined ? undefined : field.variants[key];
      return variant ? { ...normalize(variant, v), [field.tag]: key } : v;
    }
    case 'json':
      return v;
  }
}

/** A list's item: a string is trimmed (unless verbatim), then read as its shape says. */
function coerceItem(item: FieldShape, x: unknown): unknown {
  if (item.kind === 'string' && item.verbatim && typeof x === 'string') return x;
  return coerce(item, trimmed(x));
}
