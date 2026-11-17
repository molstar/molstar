/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// A language described as data (_spec/README.md "Schema"; language-description.md): its node kinds, which kinds each
// may sit under, typed params with defaults and docs. A description is plain data (no functions, no classes, a pattern
// is a string), so it can be printed, hashed, sent to an agent and turned into docs without a second description.
// Slots (`type: 'node'`) name the language they embed by id, resolved late at validation. Written from meso*'s schema
// core; `column`, `table`, `fallback` and the units hook are the delta, definitions (`def`, uses) wait.

import type { CategoryColumn, ColumnType, Node, Value } from '../tree/index.js';

export type ParamType =
  | { readonly type: 'string'; readonly enum?: readonly string[]; readonly pattern?: string }
  | {
      readonly type: 'number' | 'integer';
      readonly min?: number;
      readonly max?: number;
      /**
       * A number's unit: its bare numbers are in it, `min` and `max` too. What a quantity written in another unit
       * (`4.2nm`) means is the units hook's (`UnitsHook`): without one, a number takes no quantity.
       */
      readonly unit?: string;
      /** Whether `+inf` and `-inf` pass where it is a list's item; a single number is read by its range alone. */
      readonly infinite?: boolean;
    }
  | { readonly type: 'boolean' }
  | { readonly type: 'list'; readonly of: ParamType; readonly 'min-length'?: number; readonly 'max-length'?: number }
  | { readonly type: 'record'; readonly fields: Fields }
  /** A record whose keys are free and whose values are all of `of`; `max-keys` bounds how many. */
  | { readonly type: 'map'; readonly of: ParamType; readonly 'max-keys'?: number }
  | { readonly type: 'union'; readonly of: readonly ParamType[] }
  /** A slot: a node of another language (or this one), by the language's id. */
  | { readonly type: 'node'; readonly language: string; readonly kinds?: readonly string[] }
  | RefType
  | ColumnSpec
  | TableSpec
  | { readonly type: 'any' };

/**
 * A reference: a string naming another node of the spec, never holding it. Its value is an address word
 * (`/data/p53`, with query segments where its slot reads them) or, with `in`, a bare id resolved where `in` says.
 * Validation reads its spelling (`repr.schema.ref.address`); what it names is a resolve pass's (not built yet).
 */
export interface RefType {
  readonly type: 'ref';
  /** The kinds it may name; all when absent. */
  readonly kinds?: readonly string[];
  /** Where a bare id resolves: a kind (the nearest node of it) or an address (`/views`). */
  readonly in?: string;
}

/**
 * A typed column (format.md): one element type, read into its typed array on normalization. A float column's missing
 * value is NaN (`null` in JSON); a category's is `null` when `missing` says so; an integer column has none (`missing`
 * on one waits for its validity list).
 */
export interface ColumnSpec {
  readonly type: 'column';
  readonly of: ColumnType;
  readonly missing?: boolean;
  readonly 'min-length'?: number;
  readonly 'max-length'?: number;
}

/**
 * A table: a record of columns of equal length (each a `column` or a `list`), checked. A column is required by its own
 * `required`, as a record's field is. Extra columns are allowed and kept (a table's custom fields), each a column or a
 * list of the same length. `key` names a column whose values are unique.
 */
export interface TableSpec {
  readonly type: 'table';
  readonly columns: Fields;
  readonly key?: string;
}

/** A param: its type, its docs, and either `required` or a `default` (or neither: optional, no default). */
export type ParamSpec = ParamType & { readonly doc: string; readonly required?: boolean; readonly default?: Value };

export type Fields = Readonly<Record<string, ParamSpec>>;

/**
 * A data form's positionals: the values its `args` hold, in order, each a param spec with a name for its docs; with
 * `rest`, the last one repeats up to `max` in all.
 */
export interface DataSpec {
  readonly args: readonly ArgSpec[];
  /** The last arg repeats. */
  readonly rest?: boolean;
  /** The fewest args (default: the count of required ones). */
  readonly min?: number;
  /** The most args (default: `args.length`; with `rest`, required: work is bounded). */
  readonly max?: number;
}

export type ArgSpec = ParamSpec & { readonly name: string };

export interface KindSpec {
  readonly doc: string;
  readonly params?: Fields;
  /** The kinds this one may sit under; `[]`: only at the root (or in a slot). */
  readonly parents: readonly string[];
  /** Whether a node of this kind carries an `id`; default `optional`. */
  readonly ids?: 'required' | 'optional' | 'forbidden';
  /** A data form: its positional items are values (`args`), never an id or children. */
  readonly data?: DataSpec;
  /** A built-in kind (or the language's own word for one, a layer type) a reader without this kind uses instead. */
  readonly fallback?: string;
}

export interface Language {
  /** The id slots name it by (`seq`). */
  readonly id: string;
  readonly version: string;
  readonly doc: string;
  /** The kind of a spec's root. */
  readonly root: string;
  readonly kinds: Readonly<Record<string, KindSpec>>;
}

/**
 * The tie to a units table, a hook (tree-core.md: units stay outside the core for now). Without it, a param's `unit`
 * is documentation and a quantity string in a number slot is a type issue.
 */
export interface UnitsHook {
  /** Whether a unit id is the table's (checked by `checkLanguage`). */
  isUnit(id: string): boolean;
  /**
   * A string in a number slot of `unit` (undefined: a pure number): undefined when it is no quantity; its number in
   * the slot's unit; or why it does not fit (another dimension, a pure number).
   */
  quantity(
    text: string,
    unit: string | undefined,
  ): { readonly value: number } | { readonly message: string; readonly hint?: string } | undefined;
}

/** The kinds a language declares as data forms (`KindSpec.data`): what a reader of its specs is told. */
export function dataFormsOf(languages: Iterable<Language>): ReadonlySet<string> {
  const out = new Set<string>();
  for (const l of languages) for (const [kind, spec] of Object.entries(l.kinds)) if (spec.data) out.add(kind);
  return out;
}

// The params' static type, from the description (for builders and the language's own code)

type TypedArrayOf<C> = C extends 'int8'
  ? Int8Array
  : C extends 'uint8'
    ? Uint8Array
    : C extends 'int16'
      ? Int16Array
      : C extends 'uint16'
        ? Uint16Array
        : C extends 'int32'
          ? Int32Array
          : C extends 'uint32'
            ? Uint32Array
            : C extends 'float32'
              ? Float32Array
              : C extends 'float64'
                ? Float64Array
                : never;

/**
 * A column's value as TypeScript writes it: its typed array, or a plain array normalization types (a float column's
 * `null` is NaN; a category's `null` is missing where it says so).
 */
export type ColumnValue<T> = T extends { of: 'category'; missing: true }
  ? CategoryColumn | readonly (string | null)[]
  : T extends { of: 'category' }
    ? CategoryColumn | readonly string[]
    : T extends { of: 'float32' | 'float64' }
      ? TypedArrayOf<T['of']> | readonly (number | null)[]
      : T extends { of: infer C }
        ? TypedArrayOf<C> | readonly number[]
        : never;

/**
 * A param's value as TypeScript writes it. `Q` is what a number with a `unit` may be: `number` for the values code
 * computes with, `number | string` where a quantity is authored (the builder). Lists are read-only: a spec's values
 * are never changed in place.
 */
export type ParamValue<T, Q = number> = T extends { type: 'string'; enum: readonly (infer E)[] }
  ? E
  : T extends { type: 'string' }
    ? string
    : T extends { type: 'number' | 'integer'; unit: string }
      ? Q
      : T extends { type: 'number' | 'integer' }
        ? number
        : T extends { type: 'boolean' }
          ? boolean
          : T extends { type: 'list'; of: infer O }
            ? readonly ParamValue<O, Q>[]
            : T extends { type: 'record'; fields: infer F }
              ? ParamValues<F, Q>
              : T extends { type: 'map'; of: infer O }
                ? Readonly<Record<string, ParamValue<O, Q>>>
                : T extends { type: 'union'; of: readonly (infer O)[] }
                  ? ParamValue<O, Q>
                  : T extends { type: 'node' }
                    ? Node
                    : T extends { type: 'ref' }
                      ? string
                      : T extends { type: 'column' }
                        ? ColumnValue<T>
                        : T extends { type: 'table'; columns: infer C }
                          ? ParamValues<C, Q> & { readonly [extra: string]: Value }
                          : Value;

/**
 * A kind's params as one object type: its required ones required, the rest optional. The two mapped types are
 * flattened into one.
 */
export type ParamValues<F, Q = number> = Flat<
  { [K in keyof F as F[K] extends { required: true } ? K : never]: ParamValue<F[K], Q> } & {
    [K in keyof F as F[K] extends { required: true } ? never : K]?: ParamValue<F[K], Q>;
  }
>;

type Flat<T> = { [K in keyof T]: T[K] };

/** The params of a node of kind `K` in language `L`, as TypeScript sees them (no params: an object with no keys). */
export type ParamsOf<L extends Language, K extends keyof L['kinds']> = L['kinds'][K] extends { params: infer F }
  ? ParamValues<F>
  : Readonly<Record<string, never>>;
