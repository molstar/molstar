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
    key: Element.property(l => l.element),

    // Conformation
    x: Element.property(l => l.unit.x(l.element)),
    y: Element.property(l => l.unit.y(l.element)),
    z: Element.property(l => l.unit.z(l.element)),
    id: Element.property(l => !Unit.isAtomic(l.unit) ? notAtomic() : l.unit.conformation.atomId.value(l.element)),
    occupancy: Element.property(l => !Unit.isAtomic(l.unit) ? notAtomic() : l.unit.conformation.occupancy.value(l.element)),
    B_iso_or_equiv: Element.property(l => !Unit.isAtomic(l.unit) ? notAtomic() : l.unit.conformation.B_iso_or_equiv.value(l.element)),

    // Hierarchy
    type_symbol: Element.property(l => l.unit.hierarchy.atoms.type_symbol.value(l.element)),
    label_atom_id: Element.property(l => l.unit.hierarchy.atoms.label_atom_id.value(l.element)),
    auth_atom_id: Element.property(l => l.unit.hierarchy.atoms.auth_atom_id.value(l.element)),
    label_alt_id: Element.property(l => l.unit.hierarchy.atoms.label_alt_id.value(l.element)),
    pdbx_formal_charge: Element.property(l => l.unit.hierarchy.atoms.pdbx_formal_charge.value(l.element))
}

const residue = {
    key: Element.property(l => !Unit.isAtomic(l.unit) ? notAtomic() : l.unit.hierarchy.residueKey.value(l.unit.residueIndex[l.element])),

    group_PDB: Element.property(l => !Unit.isAtomic(l.unit) ? notAtomic() : l.unit.hierarchy.residues.group_PDB.value(l.unit.residueIndex[l.element])),
    label_comp_id: Element.property(l => !Unit.isAtomic(l.unit) ? notAtomic() : l.unit.hierarchy.residues.label_comp_id.value(l.unit.residueIndex[l.element])),
    auth_comp_id: Element.property(l => !Unit.isAtomic(l.unit) ? notAtomic() : l.unit.hierarchy.residues.auth_comp_id.value(l.unit.residueIndex[l.element])),
    label_seq_id: Element.property(l => !Unit.isAtomic(l.unit) ? notAtomic() : l.unit.hierarchy.residues.label_seq_id.value(l.unit.residueIndex[l.element])),
    auth_seq_id: Element.property(l => !Unit.isAtomic(l.unit) ? notAtomic() : l.unit.hierarchy.residues.auth_seq_id.value(l.unit.residueIndex[l.element])),
    pdbx_PDB_ins_code: Element.property(l => !Unit.isAtomic(l.unit) ? notAtomic() : l.unit.hierarchy.residues.pdbx_PDB_ins_code.value(l.unit.residueIndex[l.element]))
}

const chain = {
    key: Element.property(l => !Unit.isAtomic(l.unit) ? notAtomic() : l.unit.hierarchy.chainKey.value(l.unit.chainIndex[l.element])),

    label_asym_id: Element.property(l => !Unit.isAtomic(l.unit) ? notAtomic() : l.unit.hierarchy.chains.label_asym_id.value(l.unit.chainIndex[l.element])),
    auth_asym_id: Element.property(l => !Unit.isAtomic(l.unit) ? notAtomic() : l.unit.hierarchy.chains.auth_asym_id.value(l.unit.chainIndex[l.element])),
    label_entity_id: Element.property(l => !Unit.isAtomic(l.unit) ? notAtomic() : l.unit.hierarchy.chains.label_entity_id.value(l.unit.chainIndex[l.element]))
}

function eK(l: Element.Location) { return !Unit.isAtomic(l.unit) ? notAtomic() : l.unit.hierarchy.entityKey.value(l.unit.chainIndex[l.element]); }

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