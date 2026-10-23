/**
 * Copyright (c) 2017-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import * as ColumnHelpers from './column-helpers.js';
import { Tensor as Tensors } from '@molstar/core/math/linear-algebra/tensor';
import type { Tokens } from '@molstar/core/util/token-ranges';
import { parseInt as fastParseInt, parseFloat as fastParseFloat } from '@molstar/core/util/number-parser';
import { arrayPickIndices } from '@molstar/core/data/util/array';

interface Column<T> {
  readonly schema: Column.Schema;
  /**
   * The storage behind the values, for the fast paths of views, windows and comparisons; absent or `undefined` for a
   * computed column. It holds the values before the schema's `transform`, and a column whose `value` maps it otherwise
   * (an int column over strings) shows in `typeof value(0)`, which its readers check. May decode on its first call.
   */
  storage?(): Column.Storage | undefined;

  readonly isDefined: boolean;
  readonly rowCount: number;
  value(row: number): T;
  valueKind(row: number): Column.ValueKind;
  toArray(params?: Column.ToArrayParams<T>): ArrayLike<T>;
  areValuesEqual(rowA: number, rowB: number): boolean;
}

namespace Column {
  export type ArrayCtor<T> = { new (size: number): ArrayLike<T> };

  /**
   * What a column's values are stored in (`Column.storage`): an array of them, or a dictionary, `indices` into the
   * `values` (`value(i)` before the transform is `values[indices[i]]`; `values` may repeat a value). A reader checks
   * `typeof value(0)` against the storage first: an int column over a string field keeps the strings' storage.
   */
  export type Storage =
    | { readonly kind: 'array'; readonly array: ArrayLike<any> }
    | { readonly kind: 'dictionary'; readonly indices: ArrayLike<number>; readonly values: ArrayLike<any> };

  export type Schema<T = any> =
    | Schema.Str
    | Schema.Int
    | Schema.Float
    | Schema.Coordinate
    | Schema.Aliased<T>
    | Schema.Tensor
    | Schema.List<number | string>;

  export namespace Schema {
    // T also serves as a default value for undefined columns

    type Base<T extends string> = { valueType: T };
    export type Str = { '@type': 'str'; T: string; transform?: 'uppercase' | 'lowercase' } & Base<'str'>;
    export type Int = { '@type': 'int'; T: number } & Base<'int'>;
    export type Float = { '@type': 'float'; T: number } & Base<'float'>;
    export type Coordinate = { '@type': 'coord'; T: number } & Base<'float'>;

    export type Tensor = {
      '@type': 'tensor';
      T: Tensors.Data;
      space: Tensors.Space;
      baseType: Int | Float;
    } & Base<'tensor'>;
    export type Aliased<T> = {
      '@type': 'aliased';
      T: T;
      transform?: T extends string ? 'uppercase' | 'lowercase' : never;
    } & Base<T extends string ? 'str' : 'int'>;
    export type List<T extends number | string> = {
      '@type': 'list';
      T: T[];
      separator: string;
      itemParse: (x: string) => T;
    } & Base<'list'>;

    export const str: Str = { '@type': 'str', T: '', valueType: 'str' };
    export const ustr: Str = { '@type': 'str', T: '', valueType: 'str', transform: 'uppercase' };
    export const lstr: Str = { '@type': 'str', T: '', valueType: 'str', transform: 'lowercase' };
    export const int: Int = { '@type': 'int', T: 0, valueType: 'int' };
    export const coord: Coordinate = { '@type': 'coord', T: 0, valueType: 'float' };
    export const float: Float = { '@type': 'float', T: 0, valueType: 'float' };

    export function Str(options?: { defaultValue?: string; transform?: 'uppercase' | 'lowercase' }): Str {
      return { '@type': 'str', T: options?.defaultValue ?? '', transform: options?.transform, valueType: 'str' };
    }
    export function Int(defaultValue = 0): Int {
      return { '@type': 'int', T: defaultValue, valueType: 'int' };
    }
    export function Float(defaultValue = 0): Float {
      return { '@type': 'float', T: defaultValue, valueType: 'float' };
    }
    export function Tensor(space: Tensors.Space, baseType: Int | Float = float): Tensor {
      return { '@type': 'tensor', T: space.create(), space, valueType: 'tensor', baseType };
    }
    export function Vector(dim: number, baseType: Int | Float = float): Tensor {
      return Tensor(Tensors.Vector(dim, baseType['@type'] === 'int' ? Int32Array : Float64Array), baseType);
    }
    export function Matrix(rows: number, cols: number, baseType: Int | Float = float): Tensor {
      return Tensor(
        Tensors.ColumnMajorMatrix(rows, cols, baseType['@type'] === 'int' ? Int32Array : Float64Array),
        baseType,
      );
    }

    export function Aliased<T>(t: Str | Int): Aliased<T> {
      return t as any as Aliased<T>;
    }
    export function List<T extends number | string>(
      separator: string,
      itemParse: (x: string) => T,
      defaultValue: T[] = [],
    ): List<T> {
      return { '@type': 'list', T: defaultValue, separator, itemParse, valueType: 'list' };
    }
  }

  export interface ToArrayParams<T> {
    array?: ArrayCtor<T>;
    start?: number;
    /** Last row (exclusive) */
    end?: number;
  }

  export interface LambdaSpec<T extends Schema> {
    value: (row: number) => T['T'];
    rowCount: number;
    schema: T;
    valueKind?: (row: number) => ValueKind;
    areValuesEqual?: (rowA: number, rowB: number) => boolean;
  }

  export interface ArraySpec<T extends Schema> {
    array: ArrayLike<T['T']>;
    schema: T;
    valueKind?: (row: number) => ValueKind;
  }

  export interface DictionarySpec<T extends Schema> {
    /** indices into `values`, a typed array */
    indices: ArrayLike<number>;
    values: ArrayLike<T['T']>;
    schema: T;
    valueKind?: (row: number) => ValueKind;
  }

  export interface MapSpec<S extends Schema, T extends Schema> {
    f: (v: S['T']) => T['T'];
    schema: T;
    valueKind?: (row: number) => ValueKind;
  }

  export function is(v: any): v is Column<any> {
    return !!v && !!(v as Column<any>).schema && !!(v as Column<any>).value;
  }

  // Value kinds are accessed very of often
  // Using a const enum is an internal optimization and is defined separately to better support
  // compiling with isolatedModules flag in 3rd party use-cases of Mol*.
  export const enum ValueKinds {
    /** Defined value (= 0) */
    Present = 0,
    /** Expressed in CIF as `.` (= 1) */
    NotPresent = 1,
    /** Expressed in CIF as `?` (= 2) */
    Unknown = 2,
  }

  export const ValueKind = {
    /** Defined value (= 0) */
    Present: ValueKinds.Present,
    /** Expressed in CIF as `.` (= 1) */
    NotPresent: ValueKinds.NotPresent,
    /** Expressed in CIF as `?` (= 2) */
    Unknown: ValueKinds.Unknown,
  } as const;
  export type ValueKind = (typeof ValueKind)[keyof typeof ValueKinds];

  export function Undefined<T extends Schema>(rowCount: number, schema: T): Column<T['T']> {
    return constColumn(schema['T'], rowCount, schema, ValueKinds.NotPresent);
  }

  export function ofConst<T extends Schema>(v: T['T'], rowCount: number, type: T): Column<T['T']> {
    return constColumn(v, rowCount, type, ValueKinds.Present);
  }

  export function ofLambda<T extends Schema>(spec: LambdaSpec<T>): Column<T['T']> {
    return lambdaColumn(spec);
  }

  /** values [min, max] (i.e. include both values) */
  export function range(min: number, max: number): Column<number> {
    return ofLambda({
      value: (i) => i + min,
      rowCount: Math.max(max - min + 1, 0),
      schema: Schema.int,
    });
  }

  export function ofArray<T extends Column.Schema>(spec: Column.ArraySpec<T>): Column<T['T']> {
    return arrayColumn(spec);
  }

  /** A column of `values[indices[i]]`: the views and windows of it gather or slice the indices. */
  export function ofDictionary<T extends Column.Schema>(spec: Column.DictionarySpec<T>): Column<T['T']> {
    return dictionaryColumn(spec);
  }

  export function ofIntArray(array: ArrayLike<number>) {
    return arrayColumn({ array, schema: Schema.int });
  }

  export function ofFloatArray(array: ArrayLike<number>) {
    return arrayColumn({ array, schema: Schema.float });
  }

  export function ofStringArray(array: ArrayLike<string>) {
    return arrayColumn({ array, schema: Schema.str });
  }

  export function ofStringAliasArray<T extends string>(array: ArrayLike<T>) {
    return arrayColumn<Schema.Aliased<T>>({ array, schema: Schema.Aliased(Schema.str) });
  }

  export function ofStringListArray<T extends string>(array: ArrayLike<T[]>, separator = ',') {
    return arrayColumn<Schema.List<T>>({ array, schema: Schema.List<T>(separator, (x) => x as T) });
  }

  export function ofIntTokens(tokens: Tokens) {
    const { count, data, indices } = tokens;
    return lambdaColumn({
      value: (row: number) => fastParseInt(data, indices[2 * row], indices[2 * row + 1]) || 0,
      rowCount: count,
      schema: Schema.int,
    });
  }

  export function ofFloatTokens(tokens: Tokens) {
    const { count, data, indices } = tokens;
    return lambdaColumn({
      value: (row: number) => fastParseFloat(data, indices[2 * row], indices[2 * row + 1]) || 0,
      rowCount: count,
      schema: Schema.float,
    });
  }

  export function ofStringTokens(tokens: Tokens) {
    const { count, data, indices } = tokens;
    return lambdaColumn({
      value: (row: number) => {
        const ret = data.substring(indices[2 * row], indices[2 * row + 1]);
        if (ret === '.' || ret === '?') return '';
        return ret;
      },
      rowCount: count,
      schema: Schema.str,
    });
  }

  export function window<T>(column: Column<T>, start: number, end: number) {
    return windowColumn(column, start, end);
  }

  export function view<T>(column: Column<T>, indices: ArrayLike<number>, checkIndentity = true) {
    return columnView(column, indices, checkIndentity);
  }

  /** A map of the 1st occurence of each value. */
  export function createFirstIndexMap<T>(column: Column<T>) {
    return createFirstIndexMapOfColumn(column);
  }

  export function createIndexer<T, R extends number = number>(column: Column<T>) {
    return createIndexerOfColumn(column) as (e: T) => R;
  }

  export function mapToArray<T, S>(column: Column<T>, f: (v: T) => S, ctor?: ArrayCtor<S>): ArrayLike<S> {
    return mapToArrayImpl<T, S>(column, f, ctor || Array);
  }

  export function areEqual<T>(a: Column<T>, b: Column<T>) {
    return areColumnsEqual(a, b);
  }

  export function indicesOf<T>(c: Column<T>, test: (e: T) => boolean) {
    return columnIndicesOf(c, test);
  }

  /** Makes the column backed by an array. Useful for columns that are accessed often. */
  export function asArrayColumn<T>(c: Column<T>, array?: ArrayCtor<T>): Column<T> {
    if (c.storage?.()) return c;
    if (!c.isDefined) return Undefined(c.rowCount, c.schema) as any as Column<T>;
    return arrayColumn({ array: c.toArray({ array }), schema: c.schema, valueKind: c.valueKind });
  }

  export function copyToArray<T extends number>(c: Column<T>, array: { [k: number]: T; length: number }, offset = 0) {
    if (!c.isDefined) return;
    const s = c.storage?.();
    // the array when the column reads it as its values (as `Column.view` checks)
    const cArray = s?.kind === 'array' && c.rowCount > 0 && typeof c.value(0) === typeof s.array[0] ? s.array : void 0;
    if (cArray) {
      for (let i = 0, _i = cArray.length; i < _i; i++) array[offset + i] = cArray[i];
    } else {
      for (let i = 0, _i = c.rowCount; i < _i; i++) array[offset + i] = c.value(i);
    }
  }

  export function isIdentity<T extends number>(c: Column<T>) {
    for (let i = 0, _i = c.rowCount; i < _i; i++) {
      if (i !== c.value(i)) return false;
    }
    return true;
  }
}

export { Column };

function createFirstIndexMapOfColumn<T>(c: Column<T>): Map<T, number> {
  const map = new Map<T, number>();
  for (let i = 0, _i = c.rowCount; i < _i; i++) {
    const v = c.value(i);
    if (!map.has(v)) map.set(c.value(i), i);
  }
  return map;
}

function createIndexerOfColumn<T>(c: Column<T>): (value: T) => number {
  const map = new Map<T, number>();
  for (let i = 0, _i = c.rowCount; i < _i; i++) {
    const v = c.value(i);
    if (!map.has(v)) map.set(c.value(i), i);
  }
  return (v) => (map.has(v) ? map.get(v)! : -1);
}

function constColumn<T extends Column.Schema>(
  v: T['T'],
  rowCount: number,
  schema: T,
  valueKind: Column.ValueKind,
): Column<T['T']> {
  const value: Column<T['T']>['value'] = (row) => v;
  return {
    schema: schema,
    isDefined: valueKind === Column.ValueKinds.Present,
    rowCount,
    value,
    valueKind: (row) => valueKind,
    toArray: (params) => {
      const { array } = ColumnHelpers.createArray(rowCount, params);
      for (let i = 0, _i = array.length; i < _i; i++) array[i] = v;
      return array;
    },
    areValuesEqual: (rowA, rowB) => true,
  };
}

function lambdaColumn<T extends Column.Schema>({
  value,
  valueKind,
  areValuesEqual,
  rowCount,
  schema,
}: Column.LambdaSpec<T>): Column<T['T']> {
  return {
    schema: schema,
    isDefined: true,
    rowCount,
    value,
    valueKind: valueKind ? valueKind : (row) => Column.ValueKinds.Present,
    toArray: (params) => {
      const { array, start } = ColumnHelpers.createArray(rowCount, params);
      for (let i = 0, _i = array.length; i < _i; i++) array[i] = value(i + start);
      return array;
    },
    areValuesEqual: areValuesEqual ? areValuesEqual : (rowA, rowB) => value(rowA) === value(rowB),
  };
}

/** The values of a dictionary as a column of `schema` reads them, by `values` and schema: views and windows share them. */
const dictionaryReads = new WeakMap<object, Map<Column.Schema, ArrayLike<any>>>();

function dictionaryValues(values: ArrayLike<any>, schema: Column.Schema): ArrayLike<any> {
  if (schema.valueType !== 'str') return values;
  let bySchema = dictionaryReads.get(values);
  if (!bySchema) dictionaryReads.set(values, (bySchema = new Map()));
  let read = bySchema.get(schema);
  if (read) return read;
  const transform = (schema as Column.Schema.Str).transform;
  // as `arrayColumn` reads a str column: a non-string as text, transformed
  const asString = (v: any): string => (typeof v === 'string' ? v : `${v ?? schema.T}`);
  read =
    transform === 'lowercase'
      ? Array.from(values, (v) => asString(v).toLowerCase())
      : transform === 'uppercase'
        ? Array.from(values, (v) => asString(v).toUpperCase())
        : Array.prototype.some.call(values, (v: any) => typeof v !== 'string')
          ? Array.from(values, asString)
          : values;
  bySchema.set(schema, read);
  return read;
}

function dictionaryColumn<T extends Column.Schema>({
  indices,
  values,
  schema,
  valueKind,
}: Column.DictionarySpec<T>): Column<T['T']> {
  const rowCount = indices.length;
  const read = dictionaryValues(values, schema) as ArrayLike<T['T']>;
  const value: Column<T['T']>['value'] = (row) => read[indices[row]];
  const storage: Column.Storage = { kind: 'dictionary', indices, values };
  return {
    schema,
    storage: () => storage,
    isDefined: true,
    rowCount,
    value,
    valueKind: valueKind ? valueKind : (row) => Column.ValueKinds.Present,
    toArray: (params) => ColumnHelpers.createAndFillArray(rowCount, value, params),
    areValuesEqual: (rowA, rowB) => {
      const a = indices[rowA];
      const b = indices[rowB];
      // before the transform, as an array column compares
      return a === b || values[a] === values[b];
    },
  };
}

function arrayColumn<T extends Column.Schema>({ array, schema, valueKind }: Column.ArraySpec<T>): Column<T['T']> {
  const rowCount = array.length;
  const defaultValue = schema.T;
  const value: Column<T['T']>['value'] =
    schema.valueType === 'str'
      ? (schema as Column.Schema.Str).transform === 'lowercase'
        ? (row) => {
            const v = array[row];
            return typeof v === 'string' ? v.toLowerCase() : `${v ?? defaultValue}`.toLowerCase();
          }
        : (schema as Column.Schema.Str).transform === 'uppercase'
          ? (row) => {
              const v = array[row];
              return typeof v === 'string' ? v.toUpperCase() : `${v ?? defaultValue}`.toUpperCase();
            }
          : (row) => {
              const v = array[row];
              return typeof v === 'string' ? v : `${v ?? defaultValue}`;
            }
      : (row) => array[row];

  const isTyped = ColumnHelpers.isTypedArray(array);
  const storage: Column.Storage = { kind: 'array', array };
  return {
    schema: schema,
    storage: () => storage,
    isDefined: true,
    rowCount,
    value,
    valueKind: valueKind ? valueKind : (row) => Column.ValueKinds.Present,
    toArray:
      schema.valueType === 'str'
        ? (schema as Column.Schema.Str).transform === 'lowercase'
          ? (params) => {
              const { start, end } = ColumnHelpers.getArrayBounds(rowCount, params);
              const ret = new (params && typeof params.array !== 'undefined'
                ? params.array
                : (array as any).constructor)(end - start) as any;
              for (let i = 0, _i = end - start; i < _i; i++) {
                const v = array[start + i];
                ret[i] = typeof v === 'string' ? v.toLowerCase() : `${v ?? defaultValue}`.toLowerCase();
              }
              return ret;
            }
          : (schema as Column.Schema.Str).transform === 'uppercase'
            ? (params) => {
                const { start, end } = ColumnHelpers.getArrayBounds(rowCount, params);
                const ret = new (params && typeof params.array !== 'undefined'
                  ? params.array
                  : (array as any).constructor)(end - start) as any;
                for (let i = 0, _i = end - start; i < _i; i++) {
                  const v = array[start + i];
                  ret[i] = typeof v === 'string' ? v.toUpperCase() : `${v ?? defaultValue}`.toUpperCase();
                }
                return ret;
              }
            : (params) => {
                const { start, end } = ColumnHelpers.getArrayBounds(rowCount, params);
                const ret = new (params && typeof params.array !== 'undefined'
                  ? params.array
                  : (array as any).constructor)(end - start) as any;
                for (let i = 0, _i = end - start; i < _i; i++) {
                  const v = array[start + i];
                  ret[i] = typeof v === 'string' ? v : `${v ?? defaultValue}`;
                }
                return ret;
              }
        : isTyped
          ? (params) => ColumnHelpers.typedArrayWindow(array, params) as any as ReadonlyArray<T>
          : (params) => {
              const { start, end } = ColumnHelpers.getArrayBounds(rowCount, params);
              if (start === 0 && end === array.length) return array as ReadonlyArray<T['T']>;
              const ret = new (params && typeof params.array !== 'undefined'
                ? params.array
                : (array as any).constructor)(end - start) as any;
              for (let i = 0, _i = end - start; i < _i; i++) ret[i] = array[start + i];
              return ret;
            },
    areValuesEqual: (rowA, rowB) => array[rowA] === array[rowB],
  };
}

function windowColumn<T>(column: Column<T>, start: number, end: number): Column<T> {
  if (!column.isDefined) return Column.Undefined(end - start, column.schema);
  if (start === 0 && end === column.rowCount) return column;
  // no rows (`typedArrayWindow` would keep the last one)
  if (end <= start) return windowFull(column, start, start);
  const s = column.storage?.();
  if (s?.kind === 'array' && ColumnHelpers.isTypedArray(s.array)) return windowTyped(column, s.array, start, end);
  if (s?.kind === 'dictionary' && readsDictionary(column, s)) {
    const vk = column.valueKind;
    const indices = ColumnHelpers.isTypedArray(s.indices)
      ? ColumnHelpers.typedArrayWindow(s.indices, { start, end })
      : Array.prototype.slice.call(s.indices, start, end);
    return dictionaryColumn({ indices, values: s.values, schema: column.schema, valueKind: (row) => vk(start + row) });
  }
  return windowFull(column, start, end);
}

function windowTyped<T>(c: Column<T>, source: ArrayLike<any>, start: number, end: number): Column<T> {
  const array = ColumnHelpers.typedArrayWindow(source, { start, end });
  const vk = c.valueKind;
  const window = arrayColumn({ array, schema: c.schema, valueKind: (row) => vk(start + row) }) as any as Column<T>;
  const toArray = window.toArray;
  return {
    ...window,
    // an array of another kind from the source, which may make it otherwise than by copying (a BinaryCIF column
    // decodes straight to Float32), so a window of it reads as the whole column does
    toArray: (params) => {
      if (!params?.array || array instanceof params.array) return toArray(params);
      const { start: s, end: e } = ColumnHelpers.getArrayBounds(end - start, params);
      return c.toArray({ ...params, start: start + s, end: start + e });
    },
  };
}

function windowFull<T>(c: Column<T>, start: number, end: number): Column<T> {
  const v = c.value,
    vk = c.valueKind,
    ave = c.areValuesEqual;
  const value: Column<T>['value'] = start === 0 ? v : (row) => v(row + start);
  const rowCount = end - start;
  return {
    schema: c.schema,
    isDefined: c.isDefined,
    rowCount,
    value,
    valueKind: start === 0 ? vk : (row) => vk(row + start),
    toArray: (params) => {
      const { array, start: s } = ColumnHelpers.createArray(rowCount, params);
      for (let i = 0, _i = array.length; i < _i; i++) array[i] = v(i + s + start);
      return array;
    },
    areValuesEqual: start === 0 ? ave : (rowA, rowB) => ave(rowA + start, rowB + start),
  };
}

function isIdentity(map: ArrayLike<number>, rowCount: number) {
  if (map.length !== rowCount) return false;
  for (let i = 0, _i = map.length; i < _i; i++) {
    if (map[i] !== i) return false;
  }
  return true;
}

function columnView<T>(c: Column<T>, map: ArrayLike<number>, checkIdentity: boolean): Column<T> {
  if (c.rowCount === 0) return c;
  if (checkIdentity && isIdentity(map, c.rowCount)) return c;
  const s = c.storage?.();
  if (s?.kind === 'array' && typeof c.value(0) === typeof s.array[0]) return arrayView(c, s.array, map);
  if (s?.kind === 'dictionary' && readsDictionary(c, s)) {
    const vk = c.valueKind;
    const indices = arrayPickIndices(s.indices, map) as ArrayLike<number>;
    return dictionaryColumn({ indices, values: s.values, schema: c.schema, valueKind: (row) => vk(map[row]) });
  }
  return viewFull(c, map);
}

/** Whether the column reads its dictionary as its values (not, say, an int column over a dictionary of strings). */
function readsDictionary(c: Column<any>, s: { indices: ArrayLike<number>; values: ArrayLike<any> }) {
  return c.rowCount === 0 || typeof c.value(0) === typeof s.values[s.indices[0]];
}

function arrayView<T>(c: Column<T>, source: ArrayLike<any>, map: ArrayLike<number>): Column<T> {
  const ret = arrayPickIndices(source, map) as ArrayLike<T>;
  const vk = c.valueKind;
  return arrayColumn({ array: ret, schema: c.schema, valueKind: (row) => vk(map[row]) });
}

function viewFull<T>(c: Column<T>, map: ArrayLike<number>): Column<T> {
  const v = c.value,
    vk = c.valueKind,
    ave = c.areValuesEqual;
  const value: Column<T>['value'] = (row) => v(map[row]);
  const rowCount = map.length;
  return {
    schema: c.schema,
    isDefined: c.isDefined,
    rowCount,
    value,
    valueKind: (row) => vk(map[row]),
    toArray: (params) => {
      const { array, start } = ColumnHelpers.createArray(rowCount, params);
      for (let i = 0, _i = array.length; i < _i; i++) array[i] = v(map[i + start]);
      return array;
    },
    areValuesEqual: (rowA, rowB) => ave(map[rowA], map[rowB]),
  };
}

function mapToArrayImpl<T, S>(c: Column<T>, f: (v: T) => S, ctor: Column.ArrayCtor<S>): ArrayLike<S> {
  const ret = new ctor(c.rowCount) as any;
  for (let i = 0, _i = c.rowCount; i < _i; i++) ret[i] = f(c.value(i));
  return ret;
}

function areColumnsEqual(a: Column<any>, b: Column<any>) {
  if (a === b) return true;
  if (a.rowCount !== b.rowCount || a.isDefined !== b.isDefined || a.schema.valueType !== b.schema.valueType)
    return false;
  const sa = a.storage?.();
  const sb = b.storage?.();
  if (sa?.kind === 'array' && sb?.kind === 'array') return areArraysEqual(sa.array, sb.array, a.rowCount);
  // the views of one source: equal indices are equal values, others compare as values
  if (sa?.kind === 'dictionary' && sb?.kind === 'dictionary' && sa.values === sb.values) {
    const { indices: xs, values } = sa;
    const ys = sb.indices;
    for (let i = 0, _i = a.rowCount; i < _i; i++) {
      if (xs[i] !== ys[i] && values[xs[i]] !== values[ys[i]]) return false;
    }
    return true;
  }
  return areValuesEqual(a, b);
}

function areArraysEqual(xs: ArrayLike<any>, ys: ArrayLike<any>, rowCount: number) {
  for (let i = 0, _i = rowCount; i < _i; i++) {
    if (xs[i] !== ys[i]) return false;
  }
  return true;
}

function areValuesEqual(a: Column<any>, b: Column<any>) {
  const va = a.value,
    vb = b.value;
  for (let i = 0, _i = a.rowCount; i < _i; i++) {
    if (va(i) !== vb(i)) return false;
  }
  return true;
}

function columnIndicesOf<T>(c: Column<T>, test: (e: T) => boolean) {
  const ret = [],
    v = c.value;
  for (let i = 0, _i = c.rowCount; i < _i; i++) {
    if (test(v(i))) ret[ret.length] = i;
  }
  return ret;
}
