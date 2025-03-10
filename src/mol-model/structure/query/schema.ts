/**
 * Copyright (c) 2025 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Adam Midlik <midlik@gmail.com>
 */

import { SymmetryOperator } from '../../../mol-math/geometry';
import { MolScriptBuilder as MS } from '../../../mol-script/language/builder';
import { Expression } from '../../../mol-script/language/expression';
import { StructureElement, StructureProperties, Unit } from '../structure';

export interface StructureElementSchemaItem {
    operator_name?: string,
    label_entity_id?: string,
    label_asym_id?: string,
    auth_asym_id?: string,
    label_seq_id?: number
    auth_seq_id?: number
    label_comp_id?: string,
    auth_comp_id?: string,
    pdbx_PDB_ins_code?: string,
    beg_label_seq_id?: number
    end_label_seq_id?: number
    beg_auth_seq_id?: number
    end_auth_seq_id?: number
    label_atom_id?: string,
    auth_atom_id?: string,
    type_symbol?: string,
    atom_id?: number
    atom_index?: number
}

export type StructureElementSchema = StructureElementSchemaItem | StructureElementSchemaItem[]

// This is currently adapted from the MolViewSpec extension
// but could/should be futher optimized later, possible improvements:
//  - add optimized query that works directly on StructureElementSchema instead of converting to atoms query
//  - add more memory-efficient way to store StructureElements (e.g., struct of arrays, common prefix for multiple atoms in the same residue, etc.)

function schemaItemToExpression(item: StructureElementSchemaItem): Expression {
    const { and } = MS.core.logic;
    const { eq, gre: gte, lte } = MS.core.rel;
    const { macromolecular, ihm, core } = MS.struct.atomProperty;
    const propTests: Partial<Record<string, Expression>> = {};

    if (isDefined(item.label_entity_id)) {
        propTests['entity-test'] = eq([macromolecular.label_entity_id(), item.label_entity_id]);
    }

    const chainTests: Expression[] = [];
    if (isDefined(item.operator_name)) chainTests.push(eq([core.operatorName(), item.operator_name]));
    if (isDefined(item.label_asym_id)) chainTests.push(eq([macromolecular.label_asym_id(), item.label_asym_id]));
    if (isDefined(item.auth_asym_id)) chainTests.push(eq([macromolecular.auth_asym_id(), item.auth_asym_id]));

    if (chainTests.length === 1) {
        propTests['chain-test'] = chainTests[0];
    } else if (chainTests.length > 1) {
        propTests['chain-test'] = and(chainTests);
    }

    const residueTests: Expression[] = [];
    if (isDefined(item.label_seq_id)) {
        residueTests.push(ihm.hasSeqId({ 0: item.label_seq_id }));
    }
    if (isDefined(item.auth_seq_id)) residueTests.push(eq([macromolecular.auth_seq_id(), item.auth_seq_id]));
    if (isDefined(item.pdbx_PDB_ins_code)) residueTests.push(eq([macromolecular.pdbx_PDB_ins_code(), item.pdbx_PDB_ins_code]));

    if (isDefined(item.beg_label_seq_id) || isDefined(item.end_label_seq_id)) {
        residueTests.push(ihm.overlapsSeqIdRange({ beg: item.beg_label_seq_id, end: item.end_label_seq_id }));
    }

    if (isDefined(item.beg_auth_seq_id)) residueTests.push(gte([macromolecular.auth_seq_id(), item.beg_auth_seq_id]));
    if (isDefined(item.end_auth_seq_id)) residueTests.push(lte([macromolecular.auth_seq_id(), item.end_auth_seq_id]));
    if (residueTests.length === 1) {
        propTests['residue-test'] = residueTests[0];
    } else if (residueTests.length > 1) {
        propTests['residue-test'] = and(residueTests);
    }

    const atomTests: Expression[] = [];
    if (isDefined(item.label_comp_id)) atomTests.push(eq([macromolecular.label_comp_id(), item.label_comp_id]));
    if (isDefined(item.auth_comp_id)) atomTests.push(eq([macromolecular.auth_comp_id(), item.auth_comp_id]));
    if (isDefined(item.atom_id)) atomTests.push(eq([macromolecular.id(), item.atom_id]));
    if (isDefined(item.atom_index)) atomTests.push(eq([MS.struct.atomProperty.core.sourceIndex(), item.atom_index]));
    if (isDefined(item.label_atom_id)) atomTests.push(eq([macromolecular.label_atom_id(), item.label_atom_id]));
    if (isDefined(item.auth_atom_id)) atomTests.push(eq([macromolecular.auth_atom_id(), item.auth_atom_id]));
    if (isDefined(item.type_symbol)) atomTests.push(eq([MS.struct.atomProperty.core.elementSymbol(), item.type_symbol.toUpperCase()]));
    if (atomTests.length === 1) {
        propTests['atom-test'] = atomTests[0];
    } else if (atomTests.length > 1) {
        propTests['atom-test'] = and(atomTests);
    }

    return MS.struct.generator.atomGroups(propTests);
}

function toExpression(rows: StructureElementSchema): Expression {
    if (!Array.isArray(rows)) return schemaItemToExpression(rows as StructureElementSchemaItem);
    if (rows.length === 1) return toExpression(rows[0]);
    return unionExpression(rows.map(toExpression));
}

function unionExpression(expressions: Expression[]): Expression {
    return MS.struct.combinator.merge(expressions.map(e => MS.struct.modifier.union([e])));
}

function isDefined<T>(value: T | undefined | null): value is T {
    return value !== undefined && value !== null;
}

function forEachItem(schema: StructureElementSchema, f: (item: StructureElementSchemaItem) => void) {
    if (!Array.isArray(schema)) {
        f(schema);
    } else {
        for (const item of schema) {
            f(item);
        }
    }
}

function locationToSchemaItem(
    loc: StructureElement.Location,
    granularity: 'atom' | 'residue' | 'chain' = 'atom'
): StructureElementSchemaItem {
    // NOTE: Consider support for both auth_ and label_ prefixes

    if (Unit.isAtomic(loc.unit)) {
        let hasUniqueAtomId = true;
        let hasUniqueCompId = true;

        const { label_atom_id: atomIdCol, label_comp_id: compIdCol } = loc.unit.model.atomicHierarchy.atoms;
        const { label_seq_id: seqIdCol } = loc.unit.model.atomicHierarchy.residues;
        const { residueAtomSegments } = loc.unit.model.atomicHierarchy;
        const rI = residueAtomSegments.index[loc.element];
        const label_atom_id = StructureProperties.atom.label_atom_id(loc);
        const label_comp_id = StructureProperties.atom.label_comp_id(loc);
        for (let i = residueAtomSegments.offsets[rI], il = residueAtomSegments.offsets[rI + 1]; i < il; ++i) {
            if (i !== loc.element && atomIdCol.value(i) === label_atom_id) {
                hasUniqueAtomId = false;
            }
            if (compIdCol.value(i) !== label_comp_id) {
                hasUniqueCompId = false;
            }
        }

        const ret: StructureElementSchemaItem = {};

        if (granularity === 'residue' || granularity === 'chain' || (granularity === 'atom' && hasUniqueAtomId)) {
            ret.label_entity_id = StructureProperties.chain.label_entity_id(loc);
            ret.label_asym_id = StructureProperties.chain.label_asym_id(loc);
        }

        if (granularity === 'residue' || (granularity === 'atom' && hasUniqueAtomId)) {
            if (seqIdCol.valueKind(rI) > 0) {
                ret.label_comp_id = StructureProperties.atom.label_comp_id(loc);
            } else {
                ret.label_seq_id = seqIdCol.value(rI);
                if (!hasUniqueCompId && label_comp_id) {
                    ret.label_comp_id = label_comp_id;
                }
            }
        }

        if (granularity === 'atom') {
            if (hasUniqueAtomId) {
                ret.label_atom_id = label_atom_id;
            } else {
                ret.atom_index = StructureProperties.atom.sourceIndex(loc);
            }
        }

        if (loc.unit.conformation.operator.name !== SymmetryOperator.DefaultName && !loc.unit.conformation.operator.isIdentity) {
            ret.operator_name = loc.unit.conformation.operator.name;
        }

        return ret;
    } else {
        const ret: StructureElementSchemaItem = {
            label_entity_id: StructureProperties.chain.label_entity_id(loc),
            label_asym_id: StructureProperties.chain.label_asym_id(loc),
        };

        if (granularity === 'atom' || granularity === 'residue') {
            ret.beg_label_seq_id = StructureProperties.coarse.seq_id_begin(loc);
            ret.end_label_seq_id = StructureProperties.coarse.seq_id_end(loc);
        }

        if (loc.unit.conformation.operator.name !== SymmetryOperator.DefaultName && !loc.unit.conformation.operator.isIdentity) {
            ret.operator_name = loc.unit.conformation.operator.name;
        }

        return ret;
    }
}

export const StructureElementSchema = {
    toExpression,
    forEachItem,
    locationToSchemaItem
};