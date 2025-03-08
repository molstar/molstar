/**
 * Copyright (c) 2025 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Adam Midlik <midlik@gmail.com>
 */

import { MolScriptBuilder as MS } from '../../../mol-script/language/builder';
import { Expression } from '../../../mol-script/language/expression';

export interface StructureElementSchema {
    label_entity_id?: string,
    label_asym_id?: string,
    auth_asym_id?: string,
    label_seq_id?: number
    auth_seq_id?: number
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

export type StructureElementSchemaTable = { [K in keyof StructureElementSchema]: NonNullable<StructureElementSchema[K]>[] }

// This is currently adapted from the MolViewSpec extension
// but could/should be futher optimized later, possible improvements:
//  - add optimized query that works directly on StructureElementSchema instead of converting to atoms query
//  - add more memory-efficient way to store StructureElements (e.g., struct of arrays, common prefix for multiple atoms in the same residue, etc.)

function _structureElementSchemaToExpression(row: StructureElementSchema): Expression {
    const { and } = MS.core.logic;
    const { eq, gre: gte, lte } = MS.core.rel;
    const { macromolecular, ihm } = MS.struct.atomProperty;
    const propTests: Partial<Record<string, Expression>> = {};

    if (isDefined(row.label_entity_id)) {
        propTests['entity-test'] = eq([macromolecular.label_entity_id(), row.label_entity_id]);
    }

    const chainTests: Expression[] = [];
    if (isDefined(row.label_asym_id)) chainTests.push(eq([macromolecular.label_asym_id(), row.label_asym_id]));
    if (isDefined(row.auth_asym_id)) chainTests.push(eq([macromolecular.auth_asym_id(), row.auth_asym_id]));

    if (chainTests.length === 1) {
        propTests['chain-test'] = chainTests[0];
    } else if (chainTests.length > 1) {
        propTests['chain-test'] = and(chainTests);
    }

    const residueTests: Expression[] = [];
    if (isDefined(row.label_seq_id)) {
        residueTests.push(ihm.hasSeqId({ 0: row.label_seq_id }));
    }
    if (isDefined(row.auth_seq_id)) residueTests.push(eq([macromolecular.auth_seq_id(), row.auth_seq_id]));
    if (isDefined(row.pdbx_PDB_ins_code)) residueTests.push(eq([macromolecular.pdbx_PDB_ins_code(), row.pdbx_PDB_ins_code]));

    if (isDefined(row.beg_label_seq_id) || isDefined(row.end_label_seq_id)) {
        residueTests.push(ihm.overlapsSeqIdRange({ beg: row.beg_label_seq_id, end: row.end_label_seq_id }));
    }

    if (isDefined(row.beg_auth_seq_id)) residueTests.push(gte([macromolecular.auth_seq_id(), row.beg_auth_seq_id]));
    if (isDefined(row.end_auth_seq_id)) residueTests.push(lte([macromolecular.auth_seq_id(), row.end_auth_seq_id]));
    if (residueTests.length === 1) {
        propTests['residue-test'] = residueTests[0];
    } else if (residueTests.length > 1) {
        propTests['residue-test'] = and(residueTests);
    }

    const atomTests: Expression[] = [];
    if (isDefined(row.atom_id)) atomTests.push(eq([macromolecular.id(), row.atom_id]));
    if (isDefined(row.atom_index)) atomTests.push(eq([MS.struct.atomProperty.core.sourceIndex(), row.atom_index]));
    if (isDefined(row.label_atom_id)) atomTests.push(eq([macromolecular.label_atom_id(), row.label_atom_id]));
    if (isDefined(row.auth_atom_id)) atomTests.push(eq([macromolecular.auth_atom_id(), row.auth_atom_id]));
    if (isDefined(row.type_symbol)) atomTests.push(eq([MS.struct.atomProperty.core.elementSymbol(), row.type_symbol.toUpperCase()]));
    if (atomTests.length === 1) {
        propTests['atom-test'] = atomTests[0];
    } else if (atomTests.length > 1) {
        propTests['atom-test'] = and(atomTests);
    }

    return MS.struct.generator.atomGroups(propTests);
}

export function structureElementSchemaToExpression(rows: StructureElementSchema | readonly StructureElementSchema[]): Expression {
    if (!Array.isArray(rows)) return _structureElementSchemaToExpression(rows as StructureElementSchema);
    if (rows.length === 1) return structureElementSchemaToExpression(rows[0]);
    return unionExpression(rows.map(structureElementSchemaToExpression));
}

function unionExpression(expressions: Expression[]): Expression {
    return MS.struct.combinator.merge(expressions.map(e => MS.struct.modifier.union([e])));
}

function isDefined<T>(value: T | undefined | null): value is T {
    return value !== undefined && value !== null;
}