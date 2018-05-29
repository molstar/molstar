/**
 * Copyright (c) 2017-2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { CIFEncoder, createCIFEncoder, CIFCategory, CIFField } from 'mol-io/writer/cif'
// import { mmCIF_Schema } from 'mol-io/reader/cif/schema/mmcif'
import { Structure, Element } from '../structure'
import { Model } from '../model'
import P from '../query/properties'

interface Context {
    structure: Structure,
    model: Model
}

const atom_site_fields: CIFField<Element.Location>[] = [
    CIFField.str('group_PDB', P.residue.group_PDB),
    CIFField.int('id', P.atom.id),
    CIFField.str('type_symbol', P.atom.type_symbol as any),
    CIFField.str('label_atom_id', P.atom.label_atom_id),
    CIFField.str('label_alt_id', P.atom.label_alt_id),

    CIFField.str('label_comp_id', P.residue.label_comp_id),
    CIFField.int('label_seq_id', P.residue.label_seq_id),
    CIFField.str('pdbx_PDB_ins_code', P.residue.pdbx_PDB_ins_code),

    CIFField.str('label_asym_id', P.chain.label_asym_id),
    CIFField.str('label_entity_id', P.chain.label_entity_id),

    CIFField.float('Cartn_x', P.atom.x),
    CIFField.float('Cartn_y', P.atom.y),
    CIFField.float('Cartn_z', P.atom.z),
    CIFField.float('occupancy', P.atom.occupancy),
    CIFField.int('pdbx_formal_charge', P.atom.pdbx_formal_charge),

    CIFField.str('auth_atom_id', P.atom.auth_atom_id),
    CIFField.str('auth_comp_id', P.residue.auth_comp_id),
    CIFField.int('auth_seq_id', P.residue.auth_seq_id),
    CIFField.str('auth_asym_id', P.chain.auth_asym_id),

    CIFField.int('pdbx_PDB_model_num', P.unit.model_num),
    CIFField.str('operator_name', P.unit.operator_name)
];

function entityProvider({ model }: Context): CIFCategory {
    return CIFCategory.ofTable('entity', model.entities.data);
}

function atomSiteProvider({ structure }: Context): CIFCategory {
    return {
        data: void 0,
        name: 'atom_site',
        fields: atom_site_fields,
        rowCount: structure.elementCount,
        keys: () => structure.elementLocations()
    }
}

/** Doesn't start a data block */
export function encode_mmCIF_categories(encoder: CIFEncoder, structure: Structure) {
    const models = Structure.getModels(structure);
    if (models.length !== 1) throw 'Can\'t export stucture composed from multiple models.';
    const model = models[0];

    const ctx: Context = { structure, model };
    encoder.writeCategory(entityProvider, [ctx]);
    encoder.writeCategory(atomSiteProvider, [ctx]);
}

function to_mmCIF(name: string, structure: Structure, asBinary = false) {
    const w = createCIFEncoder({ binary: asBinary });
    w.startDataBlock(name);
    encode_mmCIF_categories(w, structure);
    return w.getData();
}

export default to_mmCIF