/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Parameter tables with MVS's names and signatures (generic/params-schema.ts there): `RequiredField`,
// `OptionalField`, `SimpleParamsSchema`, `UnionParamsSchema`. Added: a union may be open (`opts.open`: the params of
// a `<namespace>:<name>` case, whose `options` its registry entry checks).

import type { SvsRegistries } from "../registry.js";
import {
  type CheckContext,
  type FieldType,
  isNamespaced,
  isPlainObject,
  ptr,
  report,
  reportUnknown,
  show,
  typeError,
} from "./field-schema.js";

export interface RequiredField<V = unknown> {
  readonly type: FieldType<V>;
  readonly required: true;
  readonly description: string;
}
export function RequiredField<V>(type: FieldType<V>, description: string): RequiredField<V> {
  return { type, required: true, description };
}

export interface OptionalField<V = unknown> {
  readonly type: FieldType<V>;
  readonly required: false;
  /** The value a viewer uses when the param is absent; always a constant. */
  readonly default: V;
  readonly description: string;
}
export function OptionalField<V>(type: FieldType<V>, defaultValue: V, description: string): OptionalField<V> {
  return { type, required: false, default: defaultValue, description };
}

export type ParamField<V = unknown> = RequiredField<V> | OptionalField<V>;
export type Fields = { readonly [key: string]: ParamField };

export interface SimpleParamsSchema<F extends Fields = Fields> {
  readonly type: "simple";
  readonly fields: F;
}
export function SimpleParamsSchema<F extends Fields>(fields: F): SimpleParamsSchema<F> {
  return { type: "simple", fields };
}

export type Cases = { readonly [name: string]: SimpleParamsSchema };
export interface UnionParamsSchema<D extends string = string, C extends Cases = Cases> {
  readonly type: "union";
  readonly discriminator: D;
  readonly discriminatorDescription: string;
  readonly cases: C;
  /** The params of a `<namespace>:<name>` case; absent: closed. */
  readonly open: SimpleParamsSchema | null;
}
export function UnionParamsSchema<D extends string, C extends Cases>(
  discriminator: D,
  discriminatorDescription: string,
  cases: C,
  opts?: { open?: SimpleParamsSchema },
): UnionParamsSchema<D, C> {
  return { type: "union", discriminator, discriminatorDescription, cases, open: opts?.open ?? null };
}

export type ParamsSchema = SimpleParamsSchema | UnionParamsSchema;

export type ParamFieldValue<F> = F extends ParamField<infer V> ? V : never;
export type ValuesForFields<F extends Fields> = {
  -readonly [K in keyof F as F[K] extends RequiredField ? K : never]: ParamFieldValue<F[K]>;
} & { -readonly [K in keyof F as F[K] extends OptionalField ? K : never]?: ParamFieldValue<F[K]> };
/** The params a schema accepts, as MVS's `ValuesFor`. */
export type ValuesFor<P> =
  P extends SimpleParamsSchema<infer F>
    ? ValuesForFields<F>
    : P extends UnionParamsSchema<infer D, infer C>
      ? { [K in keyof C]: { [d in D]: K } & ValuesFor<C[K]> }[keyof C]
      : never;

/** The case a union's discriminator selects: a kernel case, the open case, or undefined. */
export function caseOf(schema: UnionParamsSchema, d: unknown): SimpleParamsSchema | undefined {
  if (typeof d === "string" && Object.hasOwn(schema.cases, d)) return schema.cases[d];
  return isNamespaced(d) && schema.open ? schema.open : undefined;
}

export function checkParams(schema: ParamsSchema, values: unknown, path: string, ctx: CheckContext) {
  if (!isPlainObject(values)) return typeError(ctx, path, "an object of params", values);
  if (schema.type === "simple") return checkFields(schema.fields, values, path, ctx, null);
  const d = values[schema.discriminator];
  const at = ptr(path, schema.discriminator);
  if (d === undefined) {
    return report(
      ctx,
      "svs.params.missing",
      at,
      `missing "${schema.discriminator}": ${schema.discriminatorDescription}`,
    );
  }
  const c = caseOf(schema, d);
  if (!c) {
    const known = Object.keys(schema.cases).map((k) => `"${k}"`);
    if (schema.open) known.push('"<namespace>:<name>"');
    return reportUnknown(
      ctx,
      "svs.params.type",
      at,
      `${schema.discriminator} ${show(d)} is none of ${known.join(", ")}`,
    );
  }
  checkFields(c.fields, values, path, ctx, schema.discriminator);
}

function checkFields(
  fields: Fields,
  values: Record<string, unknown>,
  path: string,
  ctx: CheckContext,
  skip: string | null,
) {
  for (const key of Object.keys(fields)) {
    const field = fields[key];
    const v = values[key];
    if (v !== undefined) field.type.check(v, ptr(path, key), ctx);
    else if (field.required)
      report(ctx, "svs.params.missing", ptr(path, key), `missing "${key}": ${field.description}`);
  }
  for (const key of Object.keys(values)) {
    if (key === skip || Object.hasOwn(fields, key)) continue;
    const names = Object.keys(fields).join(", ");
    reportUnknown(ctx, "svs.params.unknown", ptr(path, key), `unknown param "${key}"; known: ${names || "none"}`);
  }
}

const clone = <V>(v: V): V => (v !== null && typeof v === "object" ? structuredClone(v) : v);

/** Every default filled, nested values normalized (MVS's `addParamDefaults`). Unknown cases are kept as they are. */
export function normalizeParams(schema: ParamsSchema, values: Record<string, unknown>, reg: SvsRegistries) {
  const c = schema.type === "simple" ? schema : caseOf(schema, values[schema.discriminator]);
  if (!c) return { ...values };
  const out: Record<string, unknown> = { ...values };
  for (const key of Object.keys(c.fields)) {
    const field = c.fields[key];
    if (out[key] === undefined && !field.required) out[key] = clone(field.default);
    if (out[key] !== undefined && field.type.normalize) out[key] = field.type.normalize(out[key], reg);
  }
  return out;
}

/** A params schema as a field type (`object(fields)`, palettes). */
export function paramsField<V = Record<string, unknown>>(schema: ParamsSchema, name: string): FieldType<V> {
  return {
    name,
    check: (v, path, ctx) => checkParams(schema, v, path, ctx),
    normalize: (v, reg) => (isPlainObject(v) ? (normalizeParams(schema, v, reg) as V) : v),
  };
}

/** An object of these fields; without fields, any object. */
export function object<F extends Fields>(fields?: F): FieldType<ValuesFor<SimpleParamsSchema<F>>> {
  if (!fields) {
    return {
      name: "object",
      check: (v, path, ctx) => (isPlainObject(v) ? undefined : typeError(ctx, path, "object", v)),
    };
  }
  const parts = Object.entries(fields).map(([k, f]) => `${k}${f.required ? "" : "?"}: ${f.type.name}`);
  return paramsField(SimpleParamsSchema(fields), `{${parts.join(", ")}}`);
}
