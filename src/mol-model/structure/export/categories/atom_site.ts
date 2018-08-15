/**
 * Copyright (c) 2017-2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { CifWriter } from 'mol-io/writer/cif';
import { StructureElement, Structure, StructureProperties as P } from '../../structure';
import { CifExportContext } from '../mmcif';
import CifField = CifWriter.Field
import CifCategory = CifWriter.Category
import E = CifWriter.Encodings

const atom_site_fields = CifWriter.fields<StructureElement, Structure>()
    .str('group_PDB', P.residue.group_PDB)
    .index('id')
    .str('type_symbol', P.atom.type_symbol as any)
    .str('label_atom_id', P.atom.label_atom_id)

    .str('label_comp_id', P.residue.label_comp_id)
    .int('label_seq_id', P.residue.label_seq_id, {
        encoder: E.deltaRLE,
        valueKind: (k, d) => {
            const m = k.unit.model;
            return m.atomicHierarchy.residues.label_seq_id.valueKind(m.atomicHierarchy.residueAtomSegments.index[k.element]);
        }
    })
    .str('label_alt_id', P.atom.label_alt_id)
    .str('pdbx_PDB_ins_code', P.residue.pdbx_PDB_ins_code)

    .str('label_asym_id', P.chain.label_asym_id)
    .str('label_entity_id', P.chain.label_entity_id)

    .float('Cartn_x', P.atom.x, { digitCount: 3, encoder: E.fixedPoint3 })
    .float('Cartn_y', P.atom.y, { digitCount: 3, encoder: E.fixedPoint3 })
    .float('Cartn_z', P.atom.z, { digitCount: 3, encoder: E.fixedPoint3 })
    .float('occupancy', P.atom.occupancy, { digitCount: 2, encoder: E.fixedPoint2 })
    .int('pdbx_formal_charge', P.atom.pdbx_formal_charge, { encoder: E.deltaRLE })

    .str('auth_atom_id', P.atom.auth_atom_id)
    .str('auth_comp_id', P.residue.auth_comp_id)
    .int('auth_seq_id', P.residue.auth_seq_id, { encoder: E.deltaRLE })
    .str('auth_asym_id', P.chain.auth_asym_id)

    .int('pdbx_PDB_model_num', P.unit.model_num, { encoder: E.deltaRLE })
    .str('operator_name', P.unit.operator_name, {
        shouldInclude: structure => structure.units.some(u => !u.conformation.operator.isIdentity)
    })
    .getFields();

export const _atom_site: CifCategory<CifExportContext> = {
    name: 'atom_site',
    instance({ structure }: CifExportContext) {
        return {
            fields: atom_site_fields,
            data: structure,
            rowCount: structure.elementCount,
            keys: () => structure.elementLocations()
        };
    }
}

function prefixed(prefix: string, name: string) {
    return prefix ? `${prefix}_${name}` : name;
}

function mappedProp<K, D>(loc: (key: K, data: D) => StructureElement, prop: (e: StructureElement) => any) {
    return (k: K, d: D) => prop(loc(k, d));
}

export function residueIdFields<K, D>(getLocation: (key: K, data: D) => StructureElement, prefix = ''): CifField<K, D>[] {
    return CifWriter.fields<K, D>()
        .str(prefixed(prefix, `label_comp_id`), mappedProp(getLocation, P.residue.label_comp_id))
        .int(prefixed(prefix, `label_seq_id`), mappedProp(getLocation, P.residue.label_seq_id), {
            encoder: E.deltaRLE,
            valueKind: (k, d) => {
                const e = getLocation(k, d);
                const m = e.unit.model;
                return m.atomicHierarchy.residues.label_seq_id.valueKind(m.atomicHierarchy.residueAtomSegments.index[e.element]);
            }
        })
        .str(prefixed(prefix, `pdbx_PDB_ins_code`), mappedProp(getLocation, P.residue.pdbx_PDB_ins_code))

        .str(prefixed(prefix, `label_asym_id`), mappedProp(getLocation, P.chain.label_asym_id))
        .str(prefixed(prefix, `label_entity_id`), mappedProp(getLocation, P.chain.label_entity_id))

        .str(prefixed(prefix, `auth_comp_id`), mappedProp(getLocation, P.residue.auth_comp_id))
        .int(prefixed(prefix, `auth_seq_id`), mappedProp(getLocation, P.residue.auth_seq_id), { encoder: E.deltaRLE })
        .str(prefixed(prefix, `auth_asym_id`), mappedProp(getLocation, P.chain.auth_asym_id))
        .getFields();
}

export function chainIdFields<K, D>(getLocation: (key: K, data: D) => StructureElement, prefix = ''): CifField<K, D>[] {
    return CifField.build<K, D>()
        .str(prefixed(prefix, `label_asym_id`), mappedProp(getLocation, P.chain.label_asym_id))
        .str(prefixed(prefix, `label_entity_id`), mappedProp(getLocation, P.chain.label_entity_id))
        .str(prefixed(prefix, `auth_asym_id`), mappedProp(getLocation, P.chain.auth_asym_id))
        .getFields();
}