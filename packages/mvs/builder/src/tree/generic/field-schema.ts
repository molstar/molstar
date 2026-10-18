/**
 * Copyright (c) 2023-2025 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Adam Midlik <midlik@gmail.com>
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import * as iots from 'io-ts';
import { onelinerJsonString } from '@molstar/mvs-builder/json';

/** All types that can be used in tree node params.
 * Can be extended, this is just to list them all in one place and possibly catch some typing errors */
type AllowedValueTypes = string | number | boolean | null | [number, number, number] | string[] | number[] | {};

/** Type definition for a string  */
export const str = iots.string;
/** Type definition for an integer  */
export const int = iots.Integer;
/** Type definition for a float or integer number  */
export const float = iots.number;
/** Type definition for a boolean  */
export const bool = iots.boolean;
/** Type definition for a tuple, e.g. `tuple([str, int, int])`  */
export const tuple = iots.tuple;
/** Type definition for a list/array, e.g. `list(str)`  */
export const list = iots.array;
/** Type definition for a dictionary/mapping/record, e.g. `dict(str, float)` means type `{ [K in string]: number }` */
export const dict = iots.record;

/** Type definition used to create objects, e.g. `object({ name: str, age: float }, { address: str })` means type `{ name: string, age: number, address?: string }` */
export function object<P extends iots.Props, Q extends iots.Props>(
  props: P,
  optionalProps: undefined,
  name?: string,
): iots.TypeC<P>;
export function object<P extends iots.Props, Q extends iots.Props>(
  props: P,
  optionalProps: Q,
  name?: string,
): iots.IntersectionC<[iots.TypeC<P>, iots.PartialC<Q>]>;
export function object<P extends iots.Props, Q extends iots.Props>(props: P, optionalProps?: Q, name?: string) {
  if (!optionalProps) {
    return iots.type(props, name);
  }

  if (name === undefined) {
    const nameChunks = [];
    for (const key in props) {
      nameChunks.push(`${key}: ${props[key].name}`);
    }
    for (const key in optionalProps) {
      nameChunks.push(`${key}?: ${optionalProps[key].name}`);
    }
    name = `{ ${nameChunks.join(', ')} }`;
  }
  return iots.intersection([iots.type(props), iots.partial(optionalProps)], name);
}

/** Type definition used to create partial objects, e.g. `partial({ name: str, age: float })` means type `{ name?: string, age?: number }` */
export function partial<P extends iots.Props>(props: P, name?: string) {
  if (name === undefined) {
    const nameChunks = [];
    for (const key in props) {
      nameChunks.push(`${key}?: ${props[key].name}`);
    }
    name = `{ ${nameChunks.join(', ')} }`;
  }
  return iots.partial(props, name);
}

/** Type definition for union types, e.g. `union(str, int)` means string or integer */
export function union<T1 extends iots.Mixed, T2 extends iots.Mixed, TOthers extends iots.Mixed[]>(
  first: T1,
  second: T2,
  ...others: TOthers
): iots.UnionC<[T1, T2, ...TOthers]> {
  const baseTypes: iots.Mixed[] = [];
  for (const type of [first, second, ...others]) {
    if (type instanceof iots.UnionType) {
      baseTypes.push(...type.types);
    } else {
      baseTypes.push(type);
    }
  }
  return iots.union(baseTypes as any);
}

/** Type definition for nullable types, e.g. `nullable(str)` means string or `null`  */
export function nullable<V>(type: iots.Type<V>): iots.Type<V | null> {
  return union(type, iots.null);
}

/** Type definition for literal types, e.g. `literal('red', 'green', 'blue')` means 'red' or 'green' or 'blue'.
 *
 * Example usage:
 * ```
 * export type MyColor = 'red' | 'green' | 'blue';
 * export const MyColor = literal<MyColor>('red', 'green', 'blue');
 * ```
 *
 * (it looks stupid to repeat the list of values but it will result in nicer type bundle (for MolViewStories))
 */
/** Custom field types whose values are primitives (strings, numbers, booleans): their values hold nothing nested. */
export const PrimitiveFieldTypes = new WeakSet<iots.Mixed>();

export function literal<V extends string | number | boolean>(...values: V[]) {
  if (values.length === 0) {
    throw new Error(`literal type must have at least one value`);
  }
  const typeName =
    values.length === 1 ? onelinerJsonString(values[0]) : `(${values.map((v) => onelinerJsonString(v)).join(' | ')})`;
  const valueSet = new Set(values);
  const type = new iots.Type<V>(
    typeName,
    ((value: any) => valueSet.has(value)) as any,
    (value, ctx) =>
      valueSet.has(value as any)
        ? iots.success(value as V)
        : iots.failure(value, ctx, `"${value}" is not a valid value for literal type ${typeName}`),
    (value) => value,
  );
  PrimitiveFieldTypes.add(type);
  return type;
}

interface FieldBase<V extends AllowedValueTypes = any, R extends boolean = boolean> {
  /** Definition of allowed types for the field */
  type: iots.Type<V>;
  /** If `required===true`, the value must always be defined in molviewspec format (can be `null` if `type` allows it).
   * If `required===false`, the value can be ommitted (meaning that a default should be used).
   * If `type` allows `null`, the default must be `null`. */
  required: R;
  /** Description of what the field value means */
  description: string;
}

/** Schema for param field which must always be provided (has no default value) */
export interface RequiredField<V extends AllowedValueTypes = any> extends FieldBase<V> {
  required: true;
}
export function RequiredField<V extends AllowedValueTypes>(type: iots.Type<V>, description: string): RequiredField<V> {
  return { type, required: true, description };
}

/** Schema for param field which can be dropped (meaning that a default value will be used) */
export interface OptionalField<V extends AllowedValueTypes = any> extends FieldBase<V> {
  required: false;
  /** Default value for optional field.
   * If field type allows `null`, default must be `null` (this is to avoid issues in languages that do not distinguish `null` and `undefined`). */
  default: DefaultValue<V>;
}
export function OptionalField<V extends AllowedValueTypes>(
  type: iots.Type<V>,
  defaultValue: DefaultValue<V>,
  description: string,
): OptionalField<V> {
  return { type, required: false, description, default: defaultValue };
}

/** Schema for one field in params (i.e. a value in a top-level key-value pair) */
export type Field<V extends AllowedValueTypes = any> = RequiredField<V> | OptionalField<V>;

/** Type of valid default value for value type `V` (if the type allows `null`, the default must be `null`) */
type DefaultValue<V extends AllowedValueTypes> = null extends V ? null : V;

/** Type of valid value for field of type `F` (never includes `undefined`, even if field is optional) */
export type ValueFor<F extends Field | iots.Any> =
  F extends Field<infer V> ? V : F extends iots.Any ? iots.TypeOf<F> : never;

/** Return `undefined` if `value` has correct type for `field`, regardsless of if required or optional.
 * Return description of validation issues, if `value` has wrong type. */
export function fieldValidationIssues<F extends Field>(field: F, value: any): string[] | undefined {
  if (value === undefined && !field.required) return undefined; // Value undefined treated as if field not even present (unlike null)
  if (isCertainlyValid(field.type, value)) return undefined;
  const validation = field.type.decode(value);
  if (validation._tag === 'Right') {
    return undefined;
  } else {
    return reportErrors(validation.left);
  }
}

// Inlining `reportErrors` instead of `import { PathReporter } from 'io-ts/PathReporter'`;
// because it breaks Deno usage.

function reportErrors(errors: iots.Errors): string[] | undefined {
  if (errors.length === 0) return undefined;
  return errors.map(getMessage);
}

function getMessage(e: iots.ValidationError) {
  return e.message !== undefined
    ? e.message
    : `Invalid value ${stringifyError(e.value)} supplied to ${getContextPath(e.context)}`;
}

function getContextPath(context: iots.ValidationError['context']) {
  return context.map((a) => `${a.key}: ${a.type.name}`).join('/');
}

function getFunctionName(f: Function & { displayName?: string }) {
  return f.displayName || f.name || `<function ${f.length}>`;
}

function stringifyError(v: any) {
  if (typeof v === 'function') {
    return getFunctionName(v);
  }
  if (typeof v === 'number' && !isFinite(v)) {
    if (isNaN(v)) {
      return 'NaN';
    }
    return v > 0 ? 'Infinity' : '-Infinity';
  }
  return JSON.stringify(v);
}

/**
 * `type.decode(value)` would succeed, decided without io-ts's validation contexts and results (`true` only when
 * certain; `false` sends the value to `decode`, which also reports the errors). Large documents carry color maps and
 * coordinate lists of thousands of entries, where `decode` allocates per entry. Custom types are decided by `decode`,
 * once per primitive value (a document repeats a few colors many times).
 */
export function isCertainlyValid(type: iots.Mixed, value: unknown): boolean {
  return checkerOf(type)(value);
}

type Checker = (value: unknown) => boolean;
const Checkers = new WeakMap<iots.Mixed, Checker>();

/** The check of `isCertainlyValid` for `type`, made once per type (no type tests per value). */
function checkerOf(type: iots.Mixed): Checker {
  let checker = Checkers.get(type);
  if (!checker) {
    // (a placeholder while it is made, for recursive types)
    Checkers.set(type, (value) => checkerOf(type)(value));
    checker = makeChecker(type);
    Checkers.set(type, checker);
  }
  return checker;
}

function makeChecker(type: iots.Mixed): Checker {
  if (type instanceof iots.StringType) return (v) => typeof v === 'string';
  if (type instanceof iots.NumberType) return (v) => typeof v === 'number';
  if (type instanceof iots.BooleanType) return (v) => typeof v === 'boolean';
  if (type instanceof iots.NullType) return (v) => v === null;
  if (type instanceof iots.UnionType) {
    const types = type.types as iots.Mixed[];
    const checks = types.map(checkerOf);
    // (e.g. the keys of a map, `int | str`: any string passes)
    const anyString = types.some((t) => t instanceof iots.StringType);
    return (v) => {
      if (anyString && typeof v === 'string') return true;
      for (const check of checks) if (check(v)) return true;
      return false;
    };
  }
  if (type instanceof iots.RefinementType) {
    const base = checkerOf(type.type);
    const predicate = type.predicate;
    return (v) => base(v) && predicate(v);
  }
  if (type instanceof iots.ArrayType) {
    const item = checkerOf(type.type);
    return (v) => {
      if (!Array.isArray(v)) return false;
      for (let i = 0; i < v.length; i++) if (!item(v[i])) return false;
      return true;
    };
  }
  if (type instanceof iots.TupleType) {
    const items = (type.types as iots.Mixed[]).map(checkerOf);
    return (v) => {
      if (!Array.isArray(v) || v.length !== items.length) return false;
      for (let i = 0; i < items.length; i++) if (!items[i](v[i])) return false;
      return true;
    };
  }
  if (type instanceof iots.DictionaryType) {
    const domain = checkerOf(type.domain);
    const codomain = checkerOf(type.codomain);
    return (v) => {
      if (!isRecord(v)) return false;
      const keys = Object.keys(v);
      for (let i = 0; i < keys.length; i++) {
        const key = keys[i];
        if (!domain(key) || !codomain(v[key])) return false;
      }
      return true;
    };
  }
  if (type instanceof iots.InterfaceType || type instanceof iots.PartialType) {
    const partial = type instanceof iots.PartialType;
    const props = Object.entries(type.props as Record<string, iots.Mixed>).map(
      ([key, t]) => [key, checkerOf(t)] as const,
    );
    return (v) => {
      if (!isRecord(v)) return false;
      for (const [key, check] of props) {
        const x = v[key];
        if (partial && x === undefined) continue;
        if (!check(x)) return false;
      }
      return true;
    };
  }
  if (type instanceof iots.IntersectionType) {
    const checks = (type.types as iots.Mixed[]).map(checkerOf);
    return (v) => {
      for (const check of checks) if (!check(v)) return false;
      return true;
    };
  }
  // a custom type: decided by `decode`, once per primitive value; objects go to `decode` in full
  const decoded = new Map<unknown, boolean>();
  // (the last value: maps repeat a color in runs)
  let last: unknown;
  let lastValid = false;
  return (v) => {
    if (v === last) return lastValid;
    if (v !== null && typeof v === 'object') return false;
    let valid = decoded.get(v);
    if (valid === undefined) {
      if (decoded.size > 10_000) decoded.clear();
      valid = type.decode(v)._tag === 'Right';
      decoded.set(v, valid);
    }
    last = v;
    lastValid = valid;
    return valid;
  };
}

/** As io-ts's `UnknownRecord.is`. */
function isRecord(value: unknown): value is Record<string, unknown> {
  return Object.prototype.toString.call(value) === '[object Object]';
}
