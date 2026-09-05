/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { CifBlock, CifCategory, CifField } from '../../mol-io/reader/cif/data-model';
import { Mat4 } from '../../mol-math/linear-algebra';

/**
 * `_pdbx_struct_oper_list` read into one flat array of column-major 4x4 matrices.
 *
 * The generic `getMatrices`/tensor-column path allocates a `Mat4`, two tensors and twelve
 * field-name strings per row, which dominates load time for the 10^5..10^6 operators of a
 * packed particle file.
 */
export interface ParticleOperators {
    readonly count: number
    /** 16 values per operator, column-major (`j * 4 + i`). */
    readonly data: Float64Array
    /** Only set when the file deviates from serial 1-based ids; otherwise the id is `index + 1`. */
    readonly ids?: ReadonlyArray<string> | Int32Array
    readonly index?: ReadonlyMap<string, number> | ReadonlyMap<number, number>
}

/** Resolves the tensor field names once instead of per row (see `Data.tensorFieldNameGetter`). */
function getTensorFields(category: CifCategory, key: string, rank: 1 | 2): (CifField | undefined)[] {
    const { fieldNames } = category;
    const offset = (fieldNames.includes(`${key}[0]`) || fieldNames.includes(`${key}[0][0]`)) ? 0 : 1;
    const underscore = fieldNames.includes(`${key}_1`) || fieldNames.includes(`${key}_11`);

    const fields: (CifField | undefined)[] = [];
    if (rank === 1) {
        for (let i = 0; i < 3; ++i) {
            fields.push(category.getField(underscore ? `${key}_${i + offset}` : `${key}[${i + offset}]`));
        }
    } else {
        for (let i = 0; i < 3; ++i) {
            for (let j = 0; j < 3; ++j) {
                fields.push(category.getField(underscore ? `${key}_${i + offset}${j + offset}` : `${key}[${i + offset}][${j + offset}]`));
            }
        }
    }
    return fields;
}

export function readParticleOperators(block: CifBlock): ParticleOperators {
    const category = block.categories['pdbx_struct_oper_list'];
    if (!category) return { count: 0, data: new Float64Array(0) };

    const count = category.rowCount;
    const data = new Float64Array(count * 16);

    const m = getTensorFields(category, 'matrix', 2);
    const v = getTensorFields(category, 'vector', 1);

    for (let row = 0; row < count; ++row) {
        const o = row * 16;
        for (let i = 0; i < 3; ++i) {
            for (let j = 0; j < 3; ++j) {
                const f = m[i * 3 + j];
                data[o + j * 4 + i] = f ? f.float(row) : 0;
            }
            const fv = v[i];
            data[o + 12 + i] = fv ? fv.float(row) : 0;
        }
        data[o + 15] = 1;
    }

    const idField = category.getField('id');
    if (!idField) return { count, data };

    const idInfo = isSerialIdField(idField, count);
    if (idInfo.isSerial) return { count, data };

    if (idInfo.integerIds) {
        const index = new Map<number, number>();
        for (let row = 0; row < count; ++row) index.set(idInfo.integerIds[row], row);
        return { count, data, ids: idInfo.integerIds, index };
    }

    // Rare: non-serial ids (e.g. 'P', 'X0'), pay for the string table.
    const ids = new Array<string>(count);
    const index = new Map<string, number>();
    for (let row = 0; row < count; ++row) {
        const id = idField.str(row);
        ids[row] = id;
        index.set(id, row);
    }
    return { count, data, ids, index };
}

/** Perfectly serial ids require no storage; canonical integer ids with gaps use a compact array. */
function isSerialIdField(idField: CifField, count: number): { isSerial: boolean, integerIds?: Int32Array } {
    for (let row = 0; row < count; ++row) {
        const id = idField.int(row);
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
    const { data } = ops;
    const o = operatorIndex * 16;
    for (let k = 0; k < 16; ++k) out[k] = data[o + k];
    return out;
}

/** The expanded operator combinations of one `_pdbx_struct_assembly_gen` row, as operator indices. */
export interface OperatorCombinations {
    readonly count: number
    /** Number of operators per combination. */
    readonly stride: number
    /** `stride` operator indices per combination. */
    readonly indices: Int32Array
}

/**
 * Index-based counterpart of `expandOperators`, producing the same combinations in the same
 * order without allocating a `string[]` per combination.
 */
export function expandOperatorCombinations(operatorList: string[][], ops: ParticleOperators): OperatorCombinations {
    const stride = operatorList.length;
    let count = 1;
    for (let i = 0; i < stride; ++i) count *= operatorList[i].length;

    const indices = new Int32Array(count * stride);
    const current = new Int32Array(stride);
    expandInto(operatorList, ops, indices, stride, stride - 1, current, 0);
    return { count, stride, indices };
}

function expandInto(operatorList: string[][], ops: ParticleOperators, out: Int32Array, stride: number, i: number, current: Int32Array, offset: number): number {
    if (i < 0) {
        out.set(current, offset);
        return offset + stride;
    }
    const ids = operatorList[i];
    for (let j = 0, jl = ids.length; j < jl; ++j) {
        const oi = getParticleOperatorIndex(ops, ids[j]);
        if (oi < 0) throw new Error(`Operator '${ids[j]}' not found in _pdbx_struct_oper_list.`);
        current[i] = oi;
        offset = expandInto(operatorList, ops, out, stride, i - 1, current, offset);
    }
    return offset;
}

const _scratch = Mat4();

/** Multiply the operators of a combination left-to-right; the first one is outermost. */
export function composeOperatorCombination(ops: ParticleOperators, combinations: OperatorCombinations, combinationIndex: number, out: Mat4): Mat4 {
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
export function formatOperatorCombination(ids: ReadonlyArray<string> | Int32Array | undefined, combinations: OperatorCombinations, combinationIndex: number): string {
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
