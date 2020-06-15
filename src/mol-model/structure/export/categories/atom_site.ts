/**
 * Copyright (c) 2017-2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { CifWriter } from '../../../../mol-io/writer/cif';
import { StructureElement, Structure, StructureProperties as P } from '../../structure';
import { CifExportContext } from '../mmcif';
import CifField = CifWriter.Field
import CifCategory = CifWriter.Category
import E = CifWriter.Encodings

const _label_asym_id = P.chain.label_asym_id;
function atom_site_label_asym_id(e: StructureElement.Location) {
    const l = _label_asym_id(e);
    const suffix = e.unit.conformation.operator.suffix;
    if (!suffix) return l;
    return l + suffix;
}

const _auth_asym_id = P.chain.auth_asym_id;
function atom_site_auth_asym_id(e: StructureElement.Location) {
    const l = _auth_asym_id(e);
    const suffix = e.unit.conformation.operator.suffix;
    if (!suffix) return l;
    return l + suffix;
}

const atom_site_fields = CifWriter.fields<StructureElement.Location, Structure>()
    .str('group_PDB', P.residue.group_PDB)
    .index('id')
    .str('type_symbol', P.atom.type_symbol as any)
    .str('label_atom_id', P.atom.label_atom_id)

    .str('label_comp_id', P.atom.label_comp_id)
    .int('label_seq_id', P.residue.label_seq_id, {
        encoder: E.deltaRLE,
        valueKind: (k, d) => {
            const m = k.unit.model;
            return m.atomicHierarchy.residues.label_seq_id.valueKind(m.atomicHierarchy.residueAtomSegments.index[k.element]);
        }
    })
    .str('label_alt_id', P.atom.label_alt_id)
    .str('pdbx_PDB_ins_code', P.residue.pdbx_PDB_ins_code)

    .str('label_asym_id', atom_site_label_asym_id)
    .str('label_entity_id', P.chain.label_entity_id)

    .float('Cartn_x', P.atom.x, { digitCount: 3, encoder: E.fixedPoint3 })
    .float('Cartn_y', P.atom.y, { digitCount: 3, encoder: E.fixedPoint3 })
    .float('Cartn_z', P.atom.z, { digitCount: 3, encoder: E.fixedPoint3 })
    .float('occupancy', P.atom.occupancy, { digitCount: 2, encoder: E.fixedPoint2 })
    .int('pdbx_formal_charge', P.atom.pdbx_formal_charge, {
        encoder: E.deltaRLE,
        valueKind: (k, d) =>  k.unit.model.atomicHierarchy.atoms.pdbx_formal_charge.valueKind(k.element)
    })

    .str('auth_atom_id', P.atom.auth_atom_id)
    .str('auth_comp_id', P.atom.auth_comp_id)
    .int('auth_seq_id', P.residue.auth_seq_id, { encoder: E.deltaRLE })
    .str('auth_asym_id', atom_site_auth_asym_id)

    .int('pdbx_PDB_model_num', P.unit.model_num, { encoder: E.deltaRLE })
    // .str('operator_name', P.unit.operator_name, {
    //     shouldInclude: structure => structure.units.some(u => !u.conformation.operator.isIdentity)
    // })
    .getFields();

export const _atom_site: CifCategory<CifExportContext> = {
    name: 'atom_site',
    instance({ structures }: CifExportContext) {
        return {
            fields: atom_site_fields,
            source: structures.map(s => ({
                data: s,
                rowCount: s.elementCount,
                keys: () => s.elementLocations()
            }))
        };
    }
};

function prepostfixed(prefix: string | undefined, postfix: string | undefined, name: string) {
    if (prefix && postfix) return `${prefix}_${name}_${postfix}`;
    if (prefix) return `${prefix}_${name}`;
    if (postfix) return `${name}_${postfix}`;
    return name;
}

function mappedProp<K, D>(loc: (key: K, data: D) => StructureElement.Location, prop: (e: StructureElement.Location) => any) {
    return (k: K, d: D) => prop(loc(k, d));
}

function addModelNum<K, D>(fields: CifWriter.Field.Builder<K, D>, getLocation: (key: K, data: D) => StructureElement.Location, options?: IdFieldsOptions) {
    if (options && options.includeModelNum) {
        fields.int('pdbx_PDB_model_num', mappedProp(getLocation, P.unit.model_num));
    }
}

export interface IdFieldsOptions {
    prefix?: string,
    postfix?: string,
    includeModelNum?: boolean
}

export function residueIdFields<K, D>(getLocation: (key: K, data: D) => StructureElement.Location, options?: IdFieldsOptions): CifField<K, D>[] {
    const prefix = options && options.prefix, postfix = options && options.postfix;
    const ret = CifWriter.fields<K, D>()
        .str(prepostfixed(prefix, postfix, `label_comp_id`), mappedProp(getLocation, P.atom.label_comp_id))
        .int(prepostfixed(prefix, postfix, `label_seq_id`), mappedProp(getLocation, P.residue.label_seq_id), {
            encoder: E.deltaRLE,
            valueKind: (k, d) => {
                const e = getLocation(k, d);
                const m = e.unit.model;
                return m.atomicHierarchy.residues.label_seq_id.valueKind(m.atomicHierarchy.residueAtomSegments.index[e.element]);
            }
        })
        .str(prepostfixed(prefix, postfix, `pdbx_PDB_ins_code`), mappedProp(getLocation, P.residue.pdbx_PDB_ins_code))

        .str(prepostfixed(prefix, postfix, `label_asym_id`), mappedProp(getLocation, P.chain.label_asym_id))
        .str(prepostfixed(prefix, postfix, `label_entity_id`), mappedProp(getLocation, P.chain.label_entity_id))

        .str(prepostfixed(prefix, postfix, `auth_comp_id`), mappedProp(getLocation, P.atom.auth_comp_id))
        .int(prepostfixed(prefix, postfix, `auth_seq_id`), mappedProp(getLocation, P.residue.auth_seq_id), { encoder: E.deltaRLE })
        .str(prepostfixed(prefix, postfix, `auth_asym_id`), mappedProp(getLocation, P.chain.auth_asym_id));

    addModelNum(ret, getLocation, options);
    return ret.getFields();
}

export function chainIdFields<K, D>(getLocation: (key: K, data: D) => StructureElement.Location, options?: IdFieldsOptions): CifField<K, D>[] {
    const prefix = options && options.prefix, postfix = options && options.postfix;
    const ret = CifField.build<K, D>()
        .str(prepostfixed(prefix, postfix, `label_asym_id`), mappedProp(getLocation, P.chain.label_asym_id))
        .str(prepostfixed(prefix, postfix, `label_entity_id`), mappedProp(getLocation, P.chain.label_entity_id))
        .str(prepostfixed(prefix, postfix, `auth_asym_id`), mappedProp(getLocation, P.chain.auth_asym_id));

    addModelNum(ret, getLocation, options);
    return ret.getFields();
}

export function entityIdFields<K, D>(getLocation: (key: K, data: D) => StructureElement.Location, options?: IdFieldsOptions): CifField<K, D>[] {
    const prefix = options && options.prefix, postfix = options && options.postfix;
    const ret = CifField.build<K, D>()
        .str(prepostfixed(prefix, postfix, `label_entity_id`), mappedProp(getLocation, P.chain.label_entity_id));

    addModelNum(ret, getLocation, options);
    return ret.getFields();
}

export function atomIdFields<K, D>(getLocation: (key: K, data: D) => StructureElement.Location, options?: IdFieldsOptions): CifField<K, D>[] {
    const prefix = options && options.prefix, postfix = options && options.postfix;
    const ret = CifWriter.fields<K, D>()
        .str(prepostfixed(prefix, postfix, `label_atom_id`), mappedProp(getLocation, P.atom.label_atom_id))
        .str(prepostfixed(prefix, postfix, `label_comp_id`), mappedProp(getLocation, P.atom.label_comp_id))
        .int(prepostfixed(prefix, postfix, `label_seq_id`), mappedProp(getLocation, P.residue.label_seq_id), {
            encoder: E.deltaRLE,
            valueKind: (k, d) => {
                const e = getLocation(k, d);
                const m = e.unit.model;
                return m.atomicHierarchy.residues.label_seq_id.valueKind(m.atomicHierarchy.residueAtomSegments.index[e.element]);
            }
        })
        .str(prepostfixed(prefix, postfix, `label_alt_id`), mappedProp(getLocation, P.atom.label_alt_id))
        .str(prepostfixed(prefix, postfix, `pdbx_PDB_ins_code`), mappedProp(getLocation, P.residue.pdbx_PDB_ins_code))

        .str(prepostfixed(prefix, postfix, `label_asym_id`), mappedProp(getLocation, P.chain.label_asym_id))
        .str(prepostfixed(prefix, postfix, `label_entity_id`), mappedProp(getLocation, P.chain.label_entity_id))

        .str(prepostfixed(prefix, postfix, `auth_atom_id`), mappedProp(getLocation, P.atom.auth_atom_id))
        .str(prepostfixed(prefix, postfix, `auth_comp_id`), mappedProp(getLocation, P.atom.auth_comp_id))
        .int(prepostfixed(prefix, postfix, `auth_seq_id`), mappedProp(getLocation, P.residue.auth_seq_id), { encoder: E.deltaRLE })
        .str(prepostfixed(prefix, postfix, `auth_asym_id`), mappedProp(getLocation, P.chain.auth_asym_id));

    addModelNum(ret, getLocation, options);
    return ret.getFields();
}