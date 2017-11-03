/**
 * Copyright (c) 2017 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Atom } from '../structure'

const constant = {
    true: Atom.property(l => true),
    false: Atom.property(l => false),
    zero: Atom.property(l => 0)
}

const atom = {
    key: Atom.property(l => l.atom),

    // Conformation
    x: Atom.property(l => l.unit.x(l.atom)),
    y: Atom.property(l => l.unit.y(l.atom)),
    z: Atom.property(l => l.unit.z(l.atom)),
    id: Atom.property(l => l.unit.conformation.atomId.value(l.atom)),
    occupancy: Atom.property(l => l.unit.conformation.occupancy.value(l.atom)),
    B_iso_or_equiv: Atom.property(l => l.unit.conformation.B_iso_or_equiv.value(l.atom)),

    // Hierarchy
    type_symbol: Atom.property(l => l.unit.hierarchy.atoms.type_symbol.value(l.atom)),
    label_atom_id: Atom.property(l => l.unit.hierarchy.atoms.label_atom_id.value(l.atom)),
    auth_atom_id: Atom.property(l => l.unit.hierarchy.atoms.auth_atom_id.value(l.atom)),
    label_alt_id: Atom.property(l => l.unit.hierarchy.atoms.label_alt_id.value(l.atom)),
    pdbx_formal_charge: Atom.property(l => l.unit.hierarchy.atoms.pdbx_formal_charge.value(l.atom))
}

const residue = {
    key: Atom.property(l => l.unit.hierarchy.residueKey.value(l.unit.residueIndex[l.atom])),

    group_PDB: Atom.property(l => l.unit.hierarchy.residues.group_PDB.value(l.unit.residueIndex[l.atom])),
    label_comp_id: Atom.property(l => l.unit.hierarchy.residues.label_comp_id.value(l.unit.residueIndex[l.atom])),
    auth_comp_id: Atom.property(l => l.unit.hierarchy.residues.auth_comp_id.value(l.unit.residueIndex[l.atom])),
    label_seq_id: Atom.property(l => l.unit.hierarchy.residues.label_seq_id.value(l.unit.residueIndex[l.atom])),
    auth_seq_id: Atom.property(l => l.unit.hierarchy.residues.auth_seq_id.value(l.unit.residueIndex[l.atom])),
    pdbx_PDB_ins_code: Atom.property(l => l.unit.hierarchy.residues.pdbx_PDB_ins_code.value(l.unit.residueIndex[l.atom]))
}

const chain = {
    key: Atom.property(l => l.unit.hierarchy.chainKey.value(l.unit.chainIndex[l.atom])),

    label_asym_id: Atom.property(l => l.unit.hierarchy.chains.label_asym_id.value(l.unit.chainIndex[l.atom])),
    auth_asym_id: Atom.property(l => l.unit.hierarchy.chains.auth_asym_id.value(l.unit.chainIndex[l.atom])),
    label_entity_id: Atom.property(l => l.unit.hierarchy.chains.label_entity_id.value(l.unit.chainIndex[l.atom]))
}

function eK(l: Atom.Location) { return l.unit.hierarchy.entityKey.value(l.unit.chainIndex[l.atom]); }

const entity = {
    key: eK,

    id: Atom.property(l => l.unit.hierarchy.entities.id.value(eK(l))),
    type: Atom.property(l => l.unit.hierarchy.entities.type.value(eK(l))),
    src_method: Atom.property(l => l.unit.hierarchy.entities.src_method.value(eK(l))),
    pdbx_description: Atom.property(l => l.unit.hierarchy.entities.pdbx_description.value(eK(l))),
    formula_weight: Atom.property(l => l.unit.hierarchy.entities.formula_weight.value(eK(l))),
    pdbx_number_of_molecules: Atom.property(l => l.unit.hierarchy.entities.pdbx_number_of_molecules.value(eK(l))),
    details: Atom.property(l => l.unit.hierarchy.entities.details.value(eK(l))),
    pdbx_mutation: Atom.property(l => l.unit.hierarchy.entities.pdbx_mutation.value(eK(l))),
    pdbx_fragment: Atom.property(l => l.unit.hierarchy.entities.pdbx_fragment.value(eK(l))),
    pdbx_ec: Atom.property(l => l.unit.hierarchy.entities.pdbx_ec.value(eK(l)))
}

const unit = {
    operator_name: Atom.property(l => l.unit.operator.name),
    model_num: Atom.property(l => l.unit.model.modelNum)
}

const Properties = {
    constant,
    atom,
    residue,
    chain,
    entity,
    unit
}

type Properties = typeof Properties
export default Properties