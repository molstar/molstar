/**
 * Copyright (c) 2017 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Column } from 'mol-base/collections/database'
import Iterator from 'mol-base/collections/iterator'
import * as Encoder from 'mol-io/writer/cif/encoder'
import { mmCIF_Database } from 'mol-io/reader/cif/schema/mmcif'
import CIFEncoder from 'mol-io/writer/cif/encoder/text'
import { Structure, Atom, AtomSet } from '../structure'
import { Model } from '../model'
import P from '../query/properties'

interface Context {
    structure: Structure,
    model: Model
}

function str<K, D = any>(name: string, value: (k: K, d: D) => string, valueKind?: (k: K) => Column.ValueKind): Encoder.FieldDefinition<K, any> {
    return { name, type: Encoder.FieldType.Str, value, valueKind }
}

function int<K, D = any>(name: string, value: (k: K, d: D) => number, valueKind?: (k: K) => Column.ValueKind): Encoder.FieldDefinition<K, any> {
    return { name, type: Encoder.FieldType.Int, value, valueKind }
}

function float<K, D = any>(name: string, value: (k: K, d: D) => number, valueKind?: (k: K) => Column.ValueKind): Encoder.FieldDefinition<K, any> {
    return { name, type: Encoder.FieldType.Float, value, valueKind }
}

type Entity =  mmCIF_Database['entity'];

const entity: Encoder.CategoryDefinition<number, Entity> = {
    name: 'entity',
    fields: [
        str<number, Entity>('id', (i, e) => e.id.value(i)),
        str<number, Entity>('type', (i, e) => e.type.value(i)),
        str<number, Entity>('src_method', (i, e) => e.src_method.value(i)),
        str<number, Entity>('pdbx_description', (i, e) => e.pdbx_description.value(i)),
        int<number, Entity>('formula_weight', (i, e) => e.formula_weight.value(i)),
        float<number, Entity>('pdbx_number_of_molecules', (i, e) => e.pdbx_number_of_molecules.value(i)),
        str<number, Entity>('details', (i, e) => e.details.value(i)),
        str<number, Entity>('pdbx_mutation', (i, e) => e.pdbx_mutation.value(i)),
        str<number, Entity>('pdbx_fragment', (i, e) => e.pdbx_fragment.value(i)),
        str<number, Entity>('pdbx_ec', (i, e) => e.pdbx_ec.value(i)),
    ]
}

const atom_site: Encoder.CategoryDefinition<Atom.Location> = {
    name: 'atom_site',
    fields: [
        str<Atom.Location>('group_PDB', P.residue.group_PDB),
        int<Atom.Location>('id', P.atom.id),
        str<Atom.Location>('type_symbol', P.atom.type_symbol as any),
        str<Atom.Location>('label_atom_id', P.atom.label_atom_id),
        str<Atom.Location>('label_alt_id', P.atom.label_alt_id),

        str<Atom.Location>('label_comp_id', P.residue.label_comp_id),
        int<Atom.Location>('label_seq_id', P.residue.label_seq_id),
        str<Atom.Location>('pdbx_PDB_ins_code', P.residue.pdbx_PDB_ins_code),

        str<Atom.Location>('label_asym_id', P.chain.label_asym_id),
        str<Atom.Location>('label_entity_id', P.chain.label_entity_id),

        float<Atom.Location>('Cartn_x', P.atom.x),
        float<Atom.Location>('Cartn_y', P.atom.y),
        float<Atom.Location>('Cartn_z', P.atom.z),
        float<Atom.Location>('occupancy', P.atom.occupancy),
        str<Atom.Location>('pdbx_formal_charge', P.atom.pdbx_formal_charge),

        str<Atom.Location>('auth_atom_id', P.atom.auth_atom_id),
        str<Atom.Location>('auth_comp_id', P.residue.auth_comp_id),
        int<Atom.Location>('auth_seq_id', P.residue.auth_seq_id),
        str<Atom.Location>('auth_asym_id', P.chain.auth_asym_id),

        str<Atom.Location>('pdbx_operator_name', P.unit.operator_name),
    ]
};

function entityProvider({ model }: Context): Encoder.CategoryInstance {
    return {
        data: model.hierarchy.entities,
        definition: entity,
        keys: () => Iterator.Range(0, model.hierarchy.entities._rowCount - 1),
        rowCount: model.hierarchy.entities._rowCount
    }
}

function atomSiteProvider({ structure }: Context): Encoder.CategoryInstance {
    return {
        data: void 0,
        definition: atom_site,
        keys: () => Structure.atomLocationsTransient(structure),
        rowCount: AtomSet.atomCount(structure.atoms)
    }
}

function getCifString(name: string, structure: Structure) {
    const models = Structure.getModels(structure);
    if (models.length !== 1) throw 'cant export stucture composed from multiple models.';
    const model = models[0];

    const ctx: Context = { structure, model };
    const w = new CIFEncoder();

    w.startDataBlock(name);
    w.writeCategory(entityProvider, [ctx]);
    w.writeCategory(atomSiteProvider, [ctx]);
    return w.getData();
}

export default getCifString