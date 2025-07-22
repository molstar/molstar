/**
 * Copyright (c) 2023-2025 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Adam Midlik <midlik@gmail.com>
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Column } from '../../../mol-data/db';
import { ChainIndex, ElementIndex, Model, ResidueIndex, StructureElement } from '../../../mol-model/structure';
import { Expression } from '../../../mol-script/language/expression';
import { arrayExtend, filterInPlace, range } from '../../../mol-util/array';
import { AtomRanges } from './atom-ranges';
import { IndicesAndSortings, Sorting } from './indexing';
import { MVSAnnotationRow } from './schemas';
import { isAnyDefined, isDefined } from './utils';


const EmptyArray: readonly any[] = [];


/** Return atom ranges in `model` which satisfy criteria given by `row` */
export function getAtomRangesForRow(row: MVSAnnotationRow, model: Model, instanceId: string, indices: IndicesAndSortings): AtomRanges {
    if (isDefined(row.instance_id) && row.instance_id !== instanceId) return AtomRanges.empty();

    const h = model.atomicHierarchy;
    const nAtoms = h.atoms._rowCount;

    const hasAtomIds = isAnyDefined(row.atom_id, row.atom_index);
    const hasAtomFilter = isAnyDefined(row.label_atom_id, row.auth_atom_id, row.type_symbol)
        || isDefined(row.label_comp_id) && !indices.residuesByLabelCompIdIsPure
        || isDefined(row.auth_comp_id) && !indices.residuesByAuthCompIdIsPure;
    const hasResidueFilter = isAnyDefined(row.label_seq_id, row.auth_seq_id, row.pdbx_PDB_ins_code,
        row.beg_label_seq_id, row.end_label_seq_id, row.beg_auth_seq_id, row.end_auth_seq_id,
        row.label_comp_id, row.auth_comp_id);
    const hasChainFilter = isAnyDefined(row.label_asym_id, row.auth_asym_id, row.label_entity_id);

    if (hasAtomIds) {
        const theAtom = getTheAtomForRow(model, row, indices);
        return theAtom !== undefined ? AtomRanges.single(theAtom, theAtom + 1 as ElementIndex) : AtomRanges.empty();
    }

    if (!hasChainFilter && !hasResidueFilter && !hasAtomFilter) {
        return AtomRanges.single(0 as ElementIndex, nAtoms as ElementIndex);
    }

    const qualifyingChains = getQualifyingChains(model, row, indices);
    if (!hasResidueFilter && !hasAtomFilter) {
        const chainOffsets = h.chainAtomSegments.offsets;
        const ranges = AtomRanges.empty();
        for (const iChain of qualifyingChains) {
            AtomRanges.add(ranges, chainOffsets[iChain], chainOffsets[iChain + 1]);
        }
        return ranges;
    }

    const qualifyingResidues = getQualifyingResidues(model, row, indices, qualifyingChains);
    if (!hasAtomFilter) {
        const residueOffsets = h.residueAtomSegments.offsets;
        const ranges = AtomRanges.empty();
        for (const iRes of qualifyingResidues) {
            AtomRanges.add(ranges, residueOffsets[iRes], residueOffsets[iRes + 1]);
        }
        return ranges;
    }

    const qualifyingAtoms = getQualifyingAtoms(model, row, indices, qualifyingResidues);
    const ranges = AtomRanges.empty();
    for (const iAtom of qualifyingAtoms) {
        AtomRanges.add(ranges, iAtom, iAtom + 1 as ElementIndex);
    }
    return ranges;
}

/** Return atom ranges in `model` which satisfy criteria given by any of `rows` (atoms that satisfy more rows are still included only once) */
export function getAtomRangesForRows(rows: MVSAnnotationRow[], model: Model, instanceId: string, indices: IndicesAndSortings): AtomRanges {
    return AtomRanges.union(rows.map(row => getAtomRangesForRow(row, model, instanceId, indices)));
}


/** Return an array of chain indexes which satisfy criteria given by `row` */
function getQualifyingChains(model: Model, row: MVSAnnotationRow, indices: IndicesAndSortings): readonly ChainIndex[] {
    const { auth_asym_id, label_entity_id, _rowCount: nChains } = model.atomicHierarchy.chains;
    let result: readonly ChainIndex[] | undefined = undefined;
    if (isDefined(row.label_asym_id)) {
        result = indices.chainsByLabelAsymId.get(row.label_asym_id) ?? EmptyArray;
    }
    if (isDefined(row.auth_asym_id)) {
        if (result) {
            result = result.filter(i => auth_asym_id.value(i) === row.auth_asym_id);
        } else {
            result = indices.chainsByAuthAsymId.get(row.auth_asym_id) ?? EmptyArray;
        }
    }
    if (isDefined(row.label_entity_id)) {
        if (result) {
            result = result.filter(i => label_entity_id.value(i) === row.label_entity_id);
        } else {
            result = indices.chainsByLabelEntityId.get(row.label_entity_id) ?? EmptyArray;
        }
    }
    result ??= range(nChains) as ChainIndex[];
    return result;
}

/** Return an array of residue indexes which satisfy criteria given by `row` */
function getQualifyingResidues(model: Model, row: MVSAnnotationRow, indices: IndicesAndSortings, fromChains: readonly ChainIndex[]): ResidueIndex[] {
    const { label_seq_id, auth_seq_id, pdbx_PDB_ins_code } = model.atomicHierarchy.residues;
    const { label_comp_id, auth_comp_id } = model.atomicHierarchy.atoms;
    const { residueAtomSegments, chainAtomSegments } = model.atomicHierarchy;
    const { Present } = Column.ValueKind;
    const result: ResidueIndex[] = [];
    for (const iChain of fromChains) {
        let residuesHere: readonly ResidueIndex[] | undefined = undefined;
        if (isDefined(row.label_seq_id)) {
            const sorting = indices.residuesSortedByLabelSeqId.get(iChain)!;
            residuesHere = Sorting.getKeysWithValue(sorting, row.label_seq_id);
        }
        if (isDefined(row.auth_seq_id)) {
            if (residuesHere) {
                residuesHere = residuesHere.filter(i => auth_seq_id.valueKind(i) === Present && auth_seq_id.value(i) === row.auth_seq_id);
            } else {
                const sorting = indices.residuesSortedByAuthSeqId.get(iChain)!;
                residuesHere = Sorting.getKeysWithValue(sorting, row.auth_seq_id);
            }
        }
        if (isDefined(row.pdbx_PDB_ins_code)) {
            if (residuesHere) {
                residuesHere = residuesHere.filter(i => pdbx_PDB_ins_code.value(i) === row.pdbx_PDB_ins_code);
            } else {
                residuesHere = indices.residuesByInsCode.get(iChain)!.get(row.pdbx_PDB_ins_code) ?? EmptyArray;
            }
        }
        if (isDefined(row.beg_label_seq_id) || isDefined(row.end_label_seq_id)) {
            if (residuesHere) {
                if (isDefined(row.beg_label_seq_id)) {
                    residuesHere = residuesHere.filter(i => label_seq_id.valueKind(i) === Present && label_seq_id.value(i) >= row.beg_label_seq_id!);
                }
                if (isDefined(row.end_label_seq_id)) {
                    residuesHere = residuesHere.filter(i => label_seq_id.valueKind(i) === Present && label_seq_id.value(i) <= row.end_label_seq_id!);
                }
            } else {
                const sorting = indices.residuesSortedByLabelSeqId.get(iChain)!;
                residuesHere = Sorting.getKeysWithValueInRange(sorting, row.beg_label_seq_id, row.end_label_seq_id);
            }
        }
        if (isDefined(row.beg_auth_seq_id) || isDefined(row.end_auth_seq_id)) {
            if (residuesHere) {
                if (isDefined(row.beg_auth_seq_id)) {
                    residuesHere = residuesHere.filter(i => auth_seq_id.valueKind(i) === Present && auth_seq_id.value(i) >= row.beg_auth_seq_id!);
                }
                if (isDefined(row.end_auth_seq_id)) {
                    residuesHere = residuesHere.filter(i => auth_seq_id.valueKind(i) === Present && auth_seq_id.value(i) <= row.end_auth_seq_id!);
                }
            } else {
                const sorting = indices.residuesSortedByAuthSeqId.get(iChain)!;
                residuesHere = Sorting.getKeysWithValueInRange(sorting, row.beg_auth_seq_id, row.end_auth_seq_id);
            }
        }
        if (isDefined(row.label_comp_id)) {
            if (residuesHere) {
                if (indices.residuesByLabelCompIdIsPure) {
                    residuesHere = residuesHere.filter(i => label_comp_id.value(residueAtomSegments.offsets[i]) === row.label_comp_id);
                } else {
                    residuesHere = residuesHere.filter(i => {
                        for (let iAtom = residueAtomSegments.offsets[i], stop = residueAtomSegments.offsets[i + 1]; iAtom < stop; iAtom++) {
                            if (label_comp_id.value(iAtom) === row.label_comp_id) return true;
                        }
                    });
                }
            } else {
                residuesHere = indices.residuesByLabelCompId.get(iChain)!.get(row.label_comp_id) ?? EmptyArray;
            }
        }
        if (isDefined(row.auth_comp_id)) {
            if (residuesHere) {
                if (indices.residuesByAuthCompIdIsPure) {
                    residuesHere = residuesHere.filter(i => auth_comp_id.value(residueAtomSegments.offsets[i]) === row.auth_comp_id);
                } else {
                    residuesHere = residuesHere.filter(i => {
                        for (let iAtom = residueAtomSegments.offsets[i], stop = residueAtomSegments.offsets[i + 1]; iAtom < stop; iAtom++) {
                            if (auth_comp_id.value(iAtom) === row.auth_comp_id) return true;
                        }
                    });
                }
            } else {
                residuesHere = indices.residuesByAuthCompId.get(iChain)!.get(row.auth_comp_id) ?? EmptyArray;
            }
        }
        if (!residuesHere) {
            const firstResidueForChain = residueAtomSegments.index[chainAtomSegments.offsets[iChain]];
            const firstResidueAfterChain = residueAtomSegments.index[chainAtomSegments.offsets[iChain + 1] - 1] + 1;
            residuesHere = range(firstResidueForChain, firstResidueAfterChain) as ResidueIndex[];
        }
        arrayExtend(result, residuesHere);
    }
    return result;
}

/** Return an array of atom indexes which satisfy criteria given by `row` */
function getQualifyingAtoms(model: Model, row: MVSAnnotationRow, indices: IndicesAndSortings, fromResidues: readonly ResidueIndex[]): ElementIndex[] {
    const { label_atom_id, auth_atom_id, type_symbol, label_comp_id, auth_comp_id } = model.atomicHierarchy.atoms;
    const residueAtomSegments_offsets = model.atomicHierarchy.residueAtomSegments.offsets;
    const result: ElementIndex[] = [];
    for (const iRes of fromResidues) {
        const atomIdcs = range(residueAtomSegments_offsets[iRes], residueAtomSegments_offsets[iRes + 1]) as ElementIndex[];
        if (isDefined(row.label_atom_id)) {
            filterInPlace(atomIdcs, iAtom => label_atom_id.value(iAtom) === row.label_atom_id);
        }
        if (isDefined(row.auth_atom_id)) {
            filterInPlace(atomIdcs, iAtom => auth_atom_id.value(iAtom) === row.auth_atom_id);
        }
        if (isDefined(row.type_symbol)) {
            filterInPlace(atomIdcs, iAtom => type_symbol.value(iAtom) === row.type_symbol?.toUpperCase());
        }
        if (isDefined(row.label_comp_id) && !indices.residuesByLabelCompIdIsPure) {
            filterInPlace(atomIdcs, iAtom => label_comp_id.value(iAtom) === row.label_comp_id);
        }
        if (isDefined(row.auth_comp_id) && !indices.residuesByAuthCompIdIsPure) {
            filterInPlace(atomIdcs, iAtom => auth_comp_id.value(iAtom) === row.auth_comp_id);
        }
        arrayExtend(result, atomIdcs);
    }
    return result;
}

/** Return index of atom in `model` which satistfies criteria given by `row`, if any.
 * Only works when `row.atom_id` and/or `row.atom_index` is defined (otherwise use `getAtomRangesForRow`). */
function getTheAtomForRow(model: Model, row: MVSAnnotationRow, indices: IndicesAndSortings): ElementIndex | undefined {
    let iAtom: ElementIndex | undefined = undefined;
    if (!isDefined(row.atom_id) && !isDefined(row.atom_index)) throw new Error('ArgumentError: at least one of row.atom_id, row.atom_index must be defined.');
    if (isDefined(row.atom_id) && isDefined(row.atom_index)) {
        const a1 = indices.atomsById.get(row.atom_id);
        const a2 = indices.atomsByIndex.get(row.atom_index);
        if (a1 !== a2) return undefined;
        iAtom = a1;
    }
    if (isDefined(row.atom_id)) {
        iAtom = indices.atomsById.get(row.atom_id);
    }
    if (isDefined(row.atom_index)) {
        iAtom = indices.atomsByIndex.get(row.atom_index);
    }
    if (iAtom === undefined) return undefined;
    if (!atomQualifies(model, iAtom, row)) return undefined;
    return iAtom;
}

/** Return true if `iAtom`-th atom in `model` satisfies all selection criteria given by `row`. */
export function atomQualifies(model: Model, iAtom: ElementIndex, row: MVSAnnotationRow): boolean {
    const h = model.atomicHierarchy;

    const iChain = h.chainAtomSegments.index[iAtom];
    const label_asym_id = h.chains.label_asym_id.value(iChain);
    const auth_asym_id = h.chains.auth_asym_id.value(iChain);
    const label_entity_id = h.chains.label_entity_id.value(iChain);
    if (!matches(row.label_asym_id, label_asym_id)) return false;
    if (!matches(row.auth_asym_id, auth_asym_id)) return false;
    if (!matches(row.label_entity_id, label_entity_id)) return false;

    const iRes = h.residueAtomSegments.index[iAtom];
    const label_seq_id = (h.residues.label_seq_id.valueKind(iRes) === Column.ValueKind.Present) ? h.residues.label_seq_id.value(iRes) : undefined;
    const auth_seq_id = (h.residues.auth_seq_id.valueKind(iRes) === Column.ValueKind.Present) ? h.residues.auth_seq_id.value(iRes) : undefined;
    const pdbx_PDB_ins_code = h.residues.pdbx_PDB_ins_code.value(iRes);
    if (!matches(row.label_seq_id, label_seq_id)) return false;
    if (!matches(row.auth_seq_id, auth_seq_id)) return false;
    if (!matches(row.pdbx_PDB_ins_code, pdbx_PDB_ins_code)) return false;
    if (!matchesRange(row.beg_label_seq_id, row.end_label_seq_id, label_seq_id)) return false;
    if (!matchesRange(row.beg_auth_seq_id, row.end_auth_seq_id, auth_seq_id)) return false;

    const label_comp_id = h.atoms.label_comp_id.value(iAtom);
    const auth_comp_id = h.atoms.auth_comp_id.value(iAtom);
    const label_atom_id = h.atoms.label_atom_id.value(iAtom);
    const auth_atom_id = h.atoms.auth_atom_id.value(iAtom);
    const type_symbol = h.atoms.type_symbol.value(iAtom);
    const atom_id = model.atomicConformation.atomId.value(iAtom);
    const atom_index = h.atomSourceIndex.value(iAtom);
    if (!matches(row.label_comp_id, label_comp_id)) return false;
    if (!matches(row.auth_comp_id, auth_comp_id)) return false;
    if (!matches(row.label_atom_id, label_atom_id)) return false;
    if (!matches(row.auth_atom_id, auth_atom_id)) return false;
    if (!matches(row.type_symbol?.toUpperCase(), type_symbol)) return false;
    if (!matches(row.atom_id, atom_id)) return false;
    if (!matches(row.atom_index, atom_index)) return false;

    return true;
}

/** Return true if `value` equals `requiredValue` or if `requiredValue` if not defined.  */
function matches<T>(requiredValue: T | undefined | null, value: T | undefined): boolean {
    return !isDefined(requiredValue) || value === requiredValue;
}

/** Return true if `requiredMin <= value <= requiredMax`.
 * Undefined `requiredMin` behaves like negative infinity.
 * Undefined `requiredMax` behaves like positive infinity. */
function matchesRange<T>(requiredMin: T | undefined | null, requiredMax: T | undefined | null, value: T | undefined): boolean {
    if (isDefined(requiredMin) && (!isDefined(value) || value < requiredMin)) return false;
    if (isDefined(requiredMax) && (!isDefined(value) || value > requiredMax)) return false;
    return true;
}

/** Convert an annotation row into a MolScript expression */
export function rowToExpression(row: MVSAnnotationRow): Expression {
    return StructureElement.Schema.toExpression(row);
}

/** Convert multiple annotation rows into a MolScript expression.
 * (with union semantics, i.e. an atom qualifies if it qualifies for at least one of the rows) */
export function rowsToExpression(rows: readonly MVSAnnotationRow[]): Expression {
    return StructureElement.Schema.toExpression({
        items: rows as StructureElement.SchemaItem[]
    });
}

/** Data structure for an array divided into contiguous groups */
interface GroupedArray<T> {
    /** Number of groups */
    count: number,
    /** Get size of i-th group as `offsets[i+1]-offsets[i]`.
     * Get j-th element in i-th group as `grouped[offsets[i]+j]` */
    offsets: number[],
    /** Get j-th element in i-th group as `grouped[offsets[i]+j]` */
    grouped: T[],
}

/** Return row indices grouped by `row.group_id`. Rows with `row.group_id===undefined` are treated as separate groups. */
export function groupRows(rows: readonly MVSAnnotationRow[]): GroupedArray<number> {
    let counter = 0;
    const groupMap = new Map<string, number>();
    const groups: number[] = [];
    for (let i = 0; i < rows.length; i++) {
        const group_id = rows[i].group_id;
        if (!isDefined(group_id)) {
            groups.push(counter++);
        } else {
            const groupIndex = groupMap.get(group_id);
            if (groupIndex === undefined) {
                groupMap.set(group_id, counter);
                groups.push(counter);
                counter++;
            } else {
                groups.push(groupIndex);
            }
        }
    }
    const rowIndices = range(rows.length).sort((i, j) => groups[i] - groups[j]);
    const offsets: number[] = [];
    for (let i = 0; i < rows.length; i++) {
        if (i === 0 || groups[rowIndices[i]] !== groups[rowIndices[i - 1]]) offsets.push(i);
    }
    offsets.push(rowIndices.length);
    return { count: offsets.length - 1, offsets, grouped: rowIndices };
}
