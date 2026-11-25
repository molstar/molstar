/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { type CifBlock, type CifCategory, type CifField, readCifField } from '@molstar/io/reader/cif/data-model';
import { Mat4 } from '@molstar/core/math/linear-algebra';

/**
 * `_pdbx_struct_oper_list` read as its twelve columns of numbers (the binary fields' own decoded arrays), an operator's
 * matrix composed from them where asked for (`getParticleOperator`).
 *
 * The generic `getMatrices`/tensor-column path allocates a `Mat4`, two tensors and twelve
 * field-name strings per row, which dominates load time for the 10^5..10^6 operators of a
 * packed particle file. Not composed into one array of matrices either: with the columns, 224 instead of 96 bytes an
 * operator while the list is built (29 M operators: 6.6 GB)
 */
export interface ParticleOperators {
  readonly count: number;
  /**
   * The matrix's columns by its column-major index: the values at 0, 1, 2 (its first column), 4, 5, 6, 8, 9, 10, then
   * the translation at 12, 13, 14
   */
  readonly columns: readonly ArrayLike<number>[];
  /** Only set when the file deviates from serial 1-based ids; otherwise the id is `index + 1`. */
  readonly ids?: ReadonlyArray<string> | Int32Array;
  readonly index?: ReadonlyMap<string, number> | ReadonlyMap<number, number>;
}

/** Resolves the tensor field names once instead of per row (see `Data.tensorFieldNameGetter`). */
function getTensorFields(category: CifCategory, key: string, rank: 1 | 2): (CifField | undefined)[] {
  const { fieldNames } = category;
  const offset = fieldNames.includes(`${key}[0]`) || fieldNames.includes(`${key}[0][0]`) ? 0 : 1;
  const underscore = fieldNames.includes(`${key}_1`) || fieldNames.includes(`${key}_11`);

  const fields: (CifField | undefined)[] = [];
  if (rank === 1) {
    for (let i = 0; i < 3; ++i) {
      fields.push(readCifField(category, underscore ? `${key}_${i + offset}` : `${key}[${i + offset}]`));
    }
  } else {
    for (let i = 0; i < 3; ++i) {
      for (let j = 0; j < 3; ++j) {
        fields.push(
          readCifField(
            category,
            underscore ? `${key}_${i + offset}${j + offset}` : `${key}[${i + offset}][${j + offset}]`,
          ),
        );
      }
    }
  }
  return fields;
}

export function readParticleOperators(block: CifBlock): ParticleOperators {
  const category = block.categories['pdbx_struct_oper_list'];
  if (!category) return { count: 0, columns: [] };

  const count = category.rowCount;
  // whole columns (a binary field's own decoded array), not a value per accessor call: a million operators; a missing
  // one reads as 0
  let zeros: Float64Array | undefined;
  const column = (f: CifField | undefined) => f?.toFloatArray() ?? (zeros ??= new Float64Array(count));
  const [m00, m01, m02, m10, m11, m12, m20, m21, m22] = getTensorFields(category, 'matrix', 2).map(column);
  const [v0, v1, v2] = getTensorFields(category, 'vector', 1).map(column);
  // column-major: `j * 4 + i` of row i, column j
  const columns = [m00, m10, m20, m01, m11, m21, m02, m12, m22, v0, v1, v2];

  // read once, not kept decoded with the file (`readCifField`): the operators' ids and columns, millions of rows
  const idField = readCifField(category, 'id');
  if (!idField) return { count, columns };

  const idInfo = isSerialIdField(idField, count);
  if (idInfo.isSerial) return { count, columns };

  if (idInfo.integerIds) {
    const index = new Map<number, number>();
    for (let row = 0; row < count; ++row) index.set(idInfo.integerIds[row], row);
    return { count, columns, ids: idInfo.integerIds, index };
  }

  // Rare: non-serial ids (e.g. 'P', 'X0'), pay for the string table.
  const ids = new Array<string>(count);
  const index = new Map<string, number>();
  for (let row = 0; row < count; ++row) {
    const id = idField.str(row);
    ids[row] = id;
    index.set(id, row);
  }
  return { count, columns, ids, index };
}

/** Perfectly serial ids require no storage; canonical integer ids with gaps use a compact array. */
function isSerialIdField(idField: CifField, count: number): { isSerial: boolean; integerIds?: Int32Array } {
  // the whole column at once (a binary field's own array): the serial case reads every id
  const ints = idField.toIntArray();
  for (let row = 0; row < count; ++row) {
    const id = ints[row];
    if (id !== row + 1) {
      const integerIds = new Int32Array(count);
      for (let i = 0; i < row; ++i) integerIds[i] = i + 1;
      for (let i = row; i < count; ++i) {
        const integerId = idField.int(i);
        if (idField.str(i) !== `${integerId}`) return { isSerial: false };
        integerIds[i] = integerId;
      }
      return { isSerial: false, integerIds };
    }
  }
  return { isSerial: true };
}

export function getParticleOperatorId(ops: ParticleOperators, operatorIndex: number): string {
  return ops.ids ? `${ops.ids[operatorIndex]}` : `${operatorIndex + 1}`;
}

export function getParticleOperatorIndex(ops: ParticleOperators, id: string): number {
  if (ops.ids instanceof Int32Array) {
    const integerId = +id;
    if (!Number.isInteger(integerId) || `${integerId}` !== id) return -1;
    return (ops.index as ReadonlyMap<number, number>).get(integerId) ?? -1;
  }
  if (ops.index) return (ops.index as ReadonlyMap<string, number>).get(id) ?? -1;
  const i = +id - 1;
  return Number.isInteger(i) && i >= 0 && i < ops.count ? i : -1;
}

export function getParticleOperator(ops: ParticleOperators, operatorIndex: number, out: Mat4): Mat4 {
  const c = ops.columns;
  const i = operatorIndex;
  out[0] = c[0][i];
  out[1] = c[1][i];
  out[2] = c[2][i];
  out[3] = 0;
  out[4] = c[3][i];
  out[5] = c[4][i];
  out[6] = c[5][i];
  out[7] = 0;
  out[8] = c[6][i];
  out[9] = c[7][i];
  out[10] = c[8][i];
  out[11] = 0;
  out[12] = c[9][i];
  out[13] = c[10][i];
  out[14] = c[11][i];
  out[15] = 1;
  return out;
}

/** The expanded operator combinations of one `_pdbx_struct_assembly_gen` row, as operator indices. */
export interface OperatorCombinations {
  readonly count: number;
  /** Number of operators per combination. */
  readonly stride: number;
  /** `stride` operator indices per combination. */
  readonly indices: Int32Array;
}

/**
 * Index-based counterpart of `expandOperators`, producing the same combinations in the same
 * order without allocating a `string[]` per combination.
 */
export function expandOperatorCombinations(operatorList: string[][], ops: ParticleOperators): OperatorCombinations {
  return expandIndexGroups(operatorList.map((ids) => Int32Array.from(ids, (id) => getOperatorIndex(ops, id))));
}

/**
 * `expandOperatorCombinations(parseOperatorList(expression), ops)` without a string per operator: a range (`1-5000`,
 * a gen row's operators in a packed particle file) read as the indices of its ids
 */
export function parseOperatorCombinations(expression: string, ops: ParticleOperators): OperatorCombinations {
  // the groups as `parseOperatorList` reads them: '(X0)(1-5)' is [['X0'], ['1', '2', '3', '4', '5']]
  const groupRegex = /\(?([^()]+)\)?]*/g;
  const groups: Int32Array[] = [];
  let g: RegExpExecArray | null;
  while ((g = groupRegex.exec(expression))) {
    const indices: number[] = [];
    for (const e of g[1].split(',')) {
      const dashIndex = e.indexOf('-');
      if (dashIndex > 0) {
        const from = parseInt(e.substring(0, dashIndex));
        const to = parseInt(e.substring(dashIndex + 1));
        for (let id = from; id <= to; id++) indices.push(getRangeOperatorIndex(ops, id));
      } else {
        indices.push(getOperatorIndex(ops, e.trim()));
      }
    }
    groups.push(Int32Array.from(indices));
  }
  return expandIndexGroups(groups);
}

function getOperatorIndex(ops: ParticleOperators, id: string) {
  const index = getParticleOperatorIndex(ops, id);
  if (index < 0) throw new Error(`Operator '${id}' not found in _pdbx_struct_oper_list.`);
  return index;
}

/** the index of the integer id `id` of a range (`getParticleOperatorIndex` of its string) */
function getRangeOperatorIndex(ops: ParticleOperators, id: number) {
  if (!ops.ids) {
    if (id >= 1 && id <= ops.count) return id - 1;
  } else if (ops.ids instanceof Int32Array) {
    const index = (ops.index as ReadonlyMap<number, number>).get(id);
    if (index !== undefined) return index;
  } else {
    return getOperatorIndex(ops, `${id}`);
  }
  throw new Error(`Operator '${id}' not found in _pdbx_struct_oper_list.`);
}

/** the combinations of one operator of each group, as `expandOperators` orders them (the first group varies fastest) */
function expandIndexGroups(groups: Int32Array[]): OperatorCombinations {
  const stride = groups.length;
  let count = 1;
  for (let i = 0; i < stride; ++i) count *= groups[i].length;

  const indices = new Int32Array(count * stride);
  if (stride === 1) {
    indices.set(groups[0]);
    return { count, stride, indices };
  }
  const current = new Int32Array(stride);
  expandInto(groups, indices, stride, stride - 1, current, 0);
  return { count, stride, indices };
}

function expandInto(
  groups: Int32Array[],
  out: Int32Array,
  stride: number,
  i: number,
  current: Int32Array,
  offset: number,
): number {
  if (i < 0) {
    out.set(current, offset);
    return offset + stride;
  }
  const group = groups[i];
  for (let j = 0, jl = group.length; j < jl; ++j) {
    current[i] = group[j];
    offset = expandInto(groups, out, stride, i - 1, current, offset);
  }
  return offset;
}

const _scratch = Mat4();

/** Multiply the operators of a combination left-to-right; the first one is outermost. */
export function composeOperatorCombination(
  ops: ParticleOperators,
  combinations: OperatorCombinations,
  combinationIndex: number,
  out: Mat4,
): Mat4 {
  const { stride, indices } = combinations;
  const o = combinationIndex * stride;
  if (stride === 1) return getParticleOperator(ops, indices[o], out);

  Mat4.setIdentity(out);
  for (let k = 0; k < stride; ++k) {
    getParticleOperator(ops, indices[o + k], _scratch);
    Mat4.mul(out, out, _scratch);
  }
  return out;
}

/** Takes only the ids so label closures do not retain the whole matrix array. */
export function formatOperatorCombination(
  ids: ReadonlyArray<string> | Int32Array | undefined,
  combinations: OperatorCombinations,
  combinationIndex: number,
): string {
  const { stride, indices } = combinations;
  const o = combinationIndex * stride;
  if (stride === 1) return ids ? `${ids[indices[o]]}` : `${indices[o] + 1}`;
  let s = '';
  for (let k = 0; k < stride; ++k) {
    if (k > 0) s += '×';
    s += ids ? ids[indices[o + k]] : `${indices[o + k] + 1}`;
  }
  return s;
}
