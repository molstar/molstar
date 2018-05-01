/**
 * Copyright (c) 2017-2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Column } from 'mol-data/db'
import Iterator from 'mol-data/iterator'
import * as Encoder from 'mol-io/writer/cif'
// import { mmCIF_Schema } from 'mol-io/reader/cif/schema/mmcif'
import { Structure, Element, ElementSet } from '../structure'
import { Model } from '../model'
import P from '../query/properties'

interface Context {
    structure: Structure,
    model: Model
}

function str<K, D>(name: string, value: (k: K, d: D) => string, valueKind?: (k: K) => Column.ValueKind): Encoder.FieldDefinition<K, D> {
    return { name, type: Encoder.FieldType.Str, value, valueKind }
}

function int<K, D = any>(name: string, value: (k: K, d: D) => number, valueKind?: (k: K) => Column.ValueKind): Encoder.FieldDefinition<K, D> {
    return { name, type: Encoder.FieldType.Int, value, valueKind }
}

function float<K, D = any>(name: string, value: (k: K, d: D) => number, valueKind?: (k: K) => Column.ValueKind): Encoder.FieldDefinition<K, D> {
    return { name, type: Encoder.FieldType.Float, value, valueKind }
}

// function col<K, D>(name: string, c: (data: D) => Column<any>): Encoder.FieldDefinition<K, D> {
//     const kind = c['@type'].kind;
//     // TODO: matrix/vector/support
//     const type = kind === 'str' ? Encoder.FieldType.Str : kind === 'int' ? Encoder.FieldType.Int : Encoder.FieldType.Float
//     return { name, type, value, valueKind }
// }

// type Entity = Table.Columns<typeof mmCIF_Schema.entity>

// const entity: Encoder.CategoryDefinition<number, Entity> = {
//     name: 'entity',
//     fields: ofSchema(mmCIF_Schema.entity)
// }

// [
//     str('id', (i, e) => e.id.value(i)),
//     str('type', (i, e) => e.type.value(i)),
//     str('src_method', (i, e) => e.src_method.value(i)),
//     str('pdbx_description', (i, e) => e.pdbx_description.value(i)),
//     int('formula_weight', (i, e) => e.formula_weight.value(i)),
//     float('pdbx_number_of_molecules', (i, e) => e.pdbx_number_of_molecules.value(i)),
//     str('details', (i, e) => e.details.value(i)),
//     str('pdbx_mutation', (i, e) => e.pdbx_mutation.value(i)),
//     str('pdbx_fragment', (i, e) => e.pdbx_fragment.value(i)),
//     str('pdbx_ec', (i, e) => e.pdbx_ec.value(i)),
// ]

// type AtomSite = typeof mmCIF_Schema.atom_site;
// type DataSource<Key, Schema extends Table.Schema> = { [P in keyof Schema]: (key: Key) => Schema[P]['T'] }

// export const atom_site1: Partial<DataSource<Atom.Location, AtomSite>> = {
//     group_PDB: P.residue.group_PDB,
//     id: P.atom.id,
//     type_symbol: P.atom.type_symbol as any,
//     label_atom_id: P.atom.label_atom_id,
//     //...
// }

const atom_site: Encoder.CategoryDefinition<Element.Location> = {
    name: 'atom_site',
    fields: [
        str('group_PDB', P.residue.group_PDB),
        int('id', P.atom.id),
        str('type_symbol', P.atom.type_symbol as any),
        str('label_atom_id', P.atom.label_atom_id),
        str('label_alt_id', P.atom.label_alt_id),

        str('label_comp_id', P.residue.label_comp_id),
        int('label_seq_id', P.residue.label_seq_id),
        str('pdbx_PDB_ins_code', P.residue.pdbx_PDB_ins_code),

        str('label_asym_id', P.chain.label_asym_id),
        str('label_entity_id', P.chain.label_entity_id),

        float('Cartn_x', P.atom.x),
        float('Cartn_y', P.atom.y),
        float('Cartn_z', P.atom.z),
        float('occupancy', P.atom.occupancy),
        int('pdbx_formal_charge', P.atom.pdbx_formal_charge),

        str('auth_atom_id', P.atom.auth_atom_id),
        str('auth_comp_id', P.residue.auth_comp_id),
        int('auth_seq_id', P.residue.auth_seq_id),
        str('auth_asym_id', P.chain.auth_asym_id),

        int('pdbx_PDB_model_num', P.unit.model_num),
        str('pdbx_operator_name', P.unit.operator_name)
    ]
};

function entityProvider({ model }: Context): Encoder.CategoryInstance {
    return {
        data: model.entities.data,
        definition: Encoder.CategoryDefinition.ofTable('entity', model.entities.data),
        keys: () => Iterator.Range(0, model.entities.data._rowCount - 1),
        rowCount: model.entities.data._rowCount
    }
}

function atomSiteProvider({ structure }: Context): Encoder.CategoryInstance {
    return {
        data: void 0,
        definition: atom_site,
        keys: () => Structure.elementLocationsTransient(structure),
        rowCount: ElementSet.elementCount(structure.elements)
    }
}

function to_mmCIF(name: string, structure: Structure, asBinary = false) {
    const models = Structure.getModels(structure);
    if (models.length !== 1) throw 'cant export stucture composed from multiple models.';
    const model = models[0];

    const ctx: Context = { structure, model };
    const w = Encoder.create({ binary: asBinary });

    w.startDataBlock(name);
    w.writeCategory(entityProvider, [ctx]);
    w.writeCategory(atomSiteProvider, [ctx]);
    return w.getData();
}

export default to_mmCIF