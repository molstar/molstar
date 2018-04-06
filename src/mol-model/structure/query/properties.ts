/**
 * Copyright (c) 2017 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Element, Unit } from '../structure'

const constant = {
    true: Element.property(l => true),
    false: Element.property(l => false),
    zero: Element.property(l => 0)
}

function notAtomic(): never {
    throw 'Property only available for atomic models.';
}

const atom = {
    key: Element.property(l => l.atom),

    // Conformation
    x: Element.property(l => l.unit.x(l.atom)),
    y: Element.property(l => l.unit.y(l.atom)),
    z: Element.property(l => l.unit.z(l.atom)),
    id: Element.property(l => !Unit.isAtomic(l.unit) ? notAtomic() : l.unit.conformation.atomId.value(l.atom)),
    occupancy: Element.property(l => !Unit.isAtomic(l.unit) ? notAtomic() : l.unit.conformation.occupancy.value(l.atom)),
    B_iso_or_equiv: Element.property(l => !Unit.isAtomic(l.unit) ? notAtomic() : l.unit.conformation.B_iso_or_equiv.value(l.atom)),

    // Hierarchy
    type_symbol: Element.property(l => l.unit.hierarchy.atoms.type_symbol.value(l.atom)),
    label_atom_id: Element.property(l => l.unit.hierarchy.atoms.label_atom_id.value(l.atom)),
    auth_atom_id: Element.property(l => l.unit.hierarchy.atoms.auth_atom_id.value(l.atom)),
    label_alt_id: Element.property(l => l.unit.hierarchy.atoms.label_alt_id.value(l.atom)),
    pdbx_formal_charge: Element.property(l => l.unit.hierarchy.atoms.pdbx_formal_charge.value(l.atom))
}

const residue = {
    key: Element.property(l => !Unit.isAtomic(l.unit) ? notAtomic() : l.unit.hierarchy.residueKey.value(l.unit.residueIndex[l.atom])),

    group_PDB: Element.property(l => !Unit.isAtomic(l.unit) ? notAtomic() : l.unit.hierarchy.residues.group_PDB.value(l.unit.residueIndex[l.atom])),
    label_comp_id: Element.property(l => !Unit.isAtomic(l.unit) ? notAtomic() : l.unit.hierarchy.residues.label_comp_id.value(l.unit.residueIndex[l.atom])),
    auth_comp_id: Element.property(l => !Unit.isAtomic(l.unit) ? notAtomic() : l.unit.hierarchy.residues.auth_comp_id.value(l.unit.residueIndex[l.atom])),
    label_seq_id: Element.property(l => !Unit.isAtomic(l.unit) ? notAtomic() : l.unit.hierarchy.residues.label_seq_id.value(l.unit.residueIndex[l.atom])),
    auth_seq_id: Element.property(l => !Unit.isAtomic(l.unit) ? notAtomic() : l.unit.hierarchy.residues.auth_seq_id.value(l.unit.residueIndex[l.atom])),
    pdbx_PDB_ins_code: Element.property(l => !Unit.isAtomic(l.unit) ? notAtomic() : l.unit.hierarchy.residues.pdbx_PDB_ins_code.value(l.unit.residueIndex[l.atom]))
}

const chain = {
    key: Element.property(l => !Unit.isAtomic(l.unit) ? notAtomic() : l.unit.hierarchy.chainKey.value(l.unit.chainIndex[l.atom])),

    label_asym_id: Element.property(l => !Unit.isAtomic(l.unit) ? notAtomic() : l.unit.hierarchy.chains.label_asym_id.value(l.unit.chainIndex[l.atom])),
    auth_asym_id: Element.property(l => !Unit.isAtomic(l.unit) ? notAtomic() : l.unit.hierarchy.chains.auth_asym_id.value(l.unit.chainIndex[l.atom])),
    label_entity_id: Element.property(l => !Unit.isAtomic(l.unit) ? notAtomic() : l.unit.hierarchy.chains.label_entity_id.value(l.unit.chainIndex[l.atom]))
}

function eK(l: Element.Location) { return !Unit.isAtomic(l.unit) ? notAtomic() : l.unit.hierarchy.entityKey.value(l.unit.chainIndex[l.atom]); }

const entity = {
    key: eK,

    id: Element.property(l => l.unit.hierarchy.entities.id.value(eK(l))),
    type: Element.property(l => l.unit.hierarchy.entities.type.value(eK(l))),
    src_method: Element.property(l => l.unit.hierarchy.entities.src_method.value(eK(l))),
    pdbx_description: Element.property(l => l.unit.hierarchy.entities.pdbx_description.value(eK(l))),
    formula_weight: Element.property(l => l.unit.hierarchy.entities.formula_weight.value(eK(l))),
    pdbx_number_of_molecules: Element.property(l => l.unit.hierarchy.entities.pdbx_number_of_molecules.value(eK(l))),
    details: Element.property(l => l.unit.hierarchy.entities.details.value(eK(l))),
    pdbx_mutation: Element.property(l => l.unit.hierarchy.entities.pdbx_mutation.value(eK(l))),
    pdbx_fragment: Element.property(l => l.unit.hierarchy.entities.pdbx_fragment.value(eK(l))),
    pdbx_ec: Element.property(l => l.unit.hierarchy.entities.pdbx_ec.value(eK(l)))
}

const unit = {
    operator_name: Element.property(l => l.unit.operator.name),
    model_num: Element.property(l => l.unit.model.modelNum)
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