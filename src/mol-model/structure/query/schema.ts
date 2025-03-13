/**
 * Copyright (c) 2025 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Adam Midlik <midlik@gmail.com>
 */

import { MolScriptBuilder as MS } from '../../../mol-script/language/builder';
import { Expression } from '../../../mol-script/language/expression';
import { compile } from '../../../mol-script/runtime/query/base';
import { Structure, StructureElement } from '../structure';
import { QueryContext } from './context';
import { StructureSelection } from './selection';

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
    atom_id?: number,
    atom_index?: number,
}

export interface StructureElementSchemaItems {
    // The prefix is applied to each item in the list.
    // Useful for example for referencing multiple atoms in a single residue or multiple residues in a chain.
    prefix?: StructureElementSchemaItem,
    items: StructureElementSchemaItem[] | { [K in keyof StructureElementSchemaItem]: StructureElementSchemaItem[K][] }
}

export type StructureElementSchema =
    | StructureElementSchemaItem
    | StructureElementSchemaItems

function isItems(schema: StructureElementSchema): schema is StructureElementSchemaItems {
    return !!(schema as any).items;
}

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

function toExpression(schema: StructureElementSchema): Expression {
    const expressions: Expression[] = [];
    forEachItem(schema, item => expressions.push(schemaItemToExpression(item)));
    if (expressions.length === 1) return expressions[0];
    return unionExpression(expressions);
}

function unionExpression(expressions: Expression[]): Expression {
    return MS.struct.combinator.merge(expressions.map(e => MS.struct.modifier.union([e])));
}

function isDefined<T>(value: T | undefined | null): value is T {
    return value !== undefined && value !== null;
}

/**
 * Iterate over all items in a structure element schema.
 * @param schema Schema to iterate over
 * @param f Function called for each item in the schema.
 *          The value passed to the function can be mutable and should not be
 *          modified => make a copy if the value is used outside the callback.
 */
function forEachItem(schema: StructureElementSchema, f: (item: StructureElementSchemaItem) => void) {
    if (isItems(schema)) {
        if (Array.isArray(schema.items)) {
            if (schema.prefix) {
                for (const item of schema.items) {
                    f(Object.assign({}, schema.prefix, item));
                }
            } else {
                for (const item of schema.items) {
                    f(item);
                }
            }
        } else {
            const current: any = { ...schema.prefix };
            const keys = Object.keys(schema.items);
            const n = (schema.items as any)[keys[0]].length;

            for (let i = 0; i < n; i++) {
                for (const k of keys) {
                    current[k] = (schema.items as any)[k][i];
                }
                f(current);
            }
        }
    } else {
        f(schema);
    }
}

function toLoci(structure: Structure, schema: StructureElementSchema): StructureElement.Loci {
    const expr = toExpression(schema);
    const selection = compile(expr)(new QueryContext(structure));
    return StructureSelection.toLociWithSourceUnits(selection);
}

function toBundle(structure: Structure, schema: StructureElementSchema): StructureElement.Bundle {
    const loci = toLoci(structure, schema);
    return StructureElement.Bundle.fromLoci(loci);
}

export const StructureElementSchema = {
    toExpression,
    forEachItem,
    toLoci,
    toBundle,
};