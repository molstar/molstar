/**
 * Copyright (c) 2017-2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { CifWriter } from 'mol-io/writer/cif'
import { mmCIF_Schema } from 'mol-io/reader/cif/schema/mmcif'
import { Structure, Element } from '../structure'
import { Model } from '../model'
import P from '../query/properties'

interface Context {
    structure: Structure,
    model: Model
}

import CifField = CifWriter.Field
import CifCategory = CifWriter.Category

import E = CifWriter.Encodings

const atom_site_fields: CifField<Element.Location>[] = [
    CifField.str('group_PDB', P.residue.group_PDB),
    CifField.int('id', P.atom.id, { encoder: E.deltaRLE }),
    CifField.str('type_symbol', P.atom.type_symbol as any),
    CifField.str('label_atom_id', P.atom.label_atom_id),
    CifField.str('label_alt_id', P.atom.label_alt_id),

    CifField.str('label_comp_id', P.residue.label_comp_id),
    CifField.int('label_seq_id', P.residue.label_seq_id, { encoder: E.deltaRLE }),
    CifField.str('pdbx_PDB_ins_code', P.residue.pdbx_PDB_ins_code),

    CifField.str('label_asym_id', P.chain.label_asym_id),
    CifField.str('label_entity_id', P.chain.label_entity_id),

    CifField.float('Cartn_x', P.atom.x, { digitCount: 3, encoder: E.fixedPoint3 }),
    CifField.float('Cartn_y', P.atom.y, { digitCount: 3, encoder: E.fixedPoint3 }),
    CifField.float('Cartn_z', P.atom.z, { digitCount: 3, encoder: E.fixedPoint3 }),
    CifField.float('occupancy', P.atom.occupancy, { digitCount: 2, encoder: E.fixedPoint2 }),
    CifField.int('pdbx_formal_charge', P.atom.pdbx_formal_charge, { encoder: E.deltaRLE }),

    CifField.str('auth_atom_id', P.atom.auth_atom_id),
    CifField.str('auth_comp_id', P.residue.auth_comp_id),
    CifField.int('auth_seq_id', P.residue.auth_seq_id, { encoder: E.deltaRLE }),
    CifField.str('auth_asym_id', P.chain.auth_asym_id),

    CifField.int('pdbx_PDB_model_num', P.unit.model_num, { encoder: E.deltaRLE }),
    CifField.str('operator_name', P.unit.operator_name)
];

function copy_mmCif_cat(name: keyof mmCIF_Schema) {
    return ({ model }: Context) => {
        if (model.sourceData.kind !== 'mmCIF') return CifCategory.Empty;
        const table = model.sourceData.data[name];
        if (!table || !table._rowCount) return CifCategory.Empty;
        return CifCategory.ofTable(name, table);
    };
}

function _entity({ model, structure }: Context): CifCategory {
    const keys = Structure.getEntityKeys(structure);
    return CifCategory.ofTable('entity', model.entities.data, keys);
}

function _atom_site({ structure }: Context): CifCategory {
    return {
        data: void 0,
        name: 'atom_site',
        fields: atom_site_fields,
        rowCount: structure.elementCount,
        keys: () => structure.elementLocations()
    }
}

const Categories = [
    copy_mmCif_cat('entry'),
    copy_mmCif_cat('exptl'),
    copy_mmCif_cat('cell'),
    copy_mmCif_cat('symmetry'),
    _entity,
    _atom_site
]

/** Doesn't start a data block */
export function encode_mmCIF_categories(encoder: CifWriter.Encoder, structure: Structure) {
    const models = Structure.getModels(structure);
    if (models.length !== 1) throw 'Can\'t export stucture composed from multiple models.';
    const model = models[0];

    const ctx: Context[] = [{ structure, model }];

    for (const cat of Categories) {
        encoder.writeCategory(cat, ctx);
    }
}

function to_mmCIF(name: string, structure: Structure, asBinary = false) {
    const enc = CifWriter.createEncoder({ binary: asBinary });
    enc.startDataBlock(name);
    encode_mmCIF_categories(enc, structure);
    return enc.getData();
}

export default to_mmCIF