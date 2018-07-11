/**
 * Copyright (c) 2017-2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import Element from './element'
import Unit from './unit'
import { VdwRadius } from '../model/properties/atomic';

const constant = {
    true: Element.property(l => true),
    false: Element.property(l => false),
    zero: Element.property(l => 0)
}

function notAtomic(): never {
    throw 'Property only available for atomic models.';
}

function notCoarse(kind?: string): never {
    if (!!kind) throw `Property only available for coarse models (${kind}).`;
    throw `Property only available for coarse models.`;
}

// TODO: remove the type checks?

const atom = {
    key: Element.property(l => l.element),

    // Conformation
    x: Element.property(l => l.unit.conformation.x(l.element)),
    y: Element.property(l => l.unit.conformation.y(l.element)),
    z: Element.property(l => l.unit.conformation.z(l.element)),
    id: Element.property(l => !Unit.isAtomic(l.unit) ? notAtomic() : l.unit.model.atomicConformation.atomId.value(l.element)),
    occupancy: Element.property(l => !Unit.isAtomic(l.unit) ?  notAtomic() : l.unit.model.atomicConformation.occupancy.value(l.element)),
    B_iso_or_equiv: Element.property(l => !Unit.isAtomic(l.unit) ?  notAtomic() : l.unit.model.atomicConformation.B_iso_or_equiv.value(l.element)),

    // Hierarchy
    type_symbol: Element.property(l => !Unit.isAtomic(l.unit) ? notAtomic() : l.unit.model.atomicHierarchy.atoms.type_symbol.value(l.element)),
    label_atom_id: Element.property(l => !Unit.isAtomic(l.unit) ? notAtomic() : l.unit.model.atomicHierarchy.atoms.label_atom_id.value(l.element)),
    auth_atom_id: Element.property(l => !Unit.isAtomic(l.unit) ? notAtomic() : l.unit.model.atomicHierarchy.atoms.auth_atom_id.value(l.element)),
    label_alt_id: Element.property(l => !Unit.isAtomic(l.unit) ? notAtomic() : l.unit.model.atomicHierarchy.atoms.label_alt_id.value(l.element)),
    pdbx_formal_charge: Element.property(l => !Unit.isAtomic(l.unit) ? notAtomic() : l.unit.model.atomicHierarchy.atoms.pdbx_formal_charge.value(l.element)),

    // Derived
    vdw_radius: Element.property(l => !Unit.isAtomic(l.unit) ? notAtomic() : VdwRadius(l.unit.model.atomicHierarchy.atoms.type_symbol.value(l.element))),
}

const residue = {
    key: Element.property(l => !Unit.isAtomic(l.unit) ? notAtomic() : l.unit.residueIndex[l.element]),

    group_PDB: Element.property(l => !Unit.isAtomic(l.unit) ? notAtomic() : l.unit.model.atomicHierarchy.residues.group_PDB.value(l.unit.residueIndex[l.element])),
    label_comp_id: Element.property(l => !Unit.isAtomic(l.unit) ? notAtomic() : l.unit.model.atomicHierarchy.residues.label_comp_id.value(l.unit.residueIndex[l.element])),
    auth_comp_id: Element.property(l => !Unit.isAtomic(l.unit) ? notAtomic() : l.unit.model.atomicHierarchy.residues.auth_comp_id.value(l.unit.residueIndex[l.element])),
    label_seq_id: Element.property(l => !Unit.isAtomic(l.unit) ? notAtomic() : l.unit.model.atomicHierarchy.residues.label_seq_id.value(l.unit.residueIndex[l.element])),
    auth_seq_id: Element.property(l => !Unit.isAtomic(l.unit) ? notAtomic() : l.unit.model.atomicHierarchy.residues.auth_seq_id.value(l.unit.residueIndex[l.element])),
    pdbx_PDB_ins_code: Element.property(l => !Unit.isAtomic(l.unit) ? notAtomic() : l.unit.model.atomicHierarchy.residues.pdbx_PDB_ins_code.value(l.unit.residueIndex[l.element])),

    // Properties
    secondary_structure_type: Element.property(l => !Unit.isAtomic(l.unit) ? notAtomic() : l.unit.model.properties.secondaryStructure.type[l.unit.residueIndex[l.element]]),
    secondary_structure_key: Element.property(l => !Unit.isAtomic(l.unit) ? notAtomic() : l.unit.model.properties.secondaryStructure.key[l.unit.residueIndex[l.element]]),
}

const chain = {
    key: Element.property(l => !Unit.isAtomic(l.unit) ? notAtomic() : l.unit.chainIndex[l.element]),

    label_asym_id: Element.property(l => !Unit.isAtomic(l.unit) ? notAtomic() : l.unit.model.atomicHierarchy.chains.label_asym_id.value(l.unit.chainIndex[l.element])),
    auth_asym_id: Element.property(l => !Unit.isAtomic(l.unit) ? notAtomic() : l.unit.model.atomicHierarchy.chains.auth_asym_id.value(l.unit.chainIndex[l.element])),
    label_entity_id: Element.property(l => !Unit.isAtomic(l.unit) ? notAtomic() : l.unit.model.atomicHierarchy.chains.label_entity_id.value(l.unit.chainIndex[l.element]))
}

const coarse = {
    key: atom.key,
    entityKey: Element.property(l => !Unit.isCoarse(l.unit) ? notCoarse() : l.unit.coarseElements.entityKey[l.element]),

    x: atom.x,
    y: atom.y,
    z: atom.z,

    asym_id: Element.property(l => !Unit.isCoarse(l.unit) ? notCoarse() : l.unit.coarseElements.asym_id.value(l.element)),
    seq_id_begin: Element.property(l => !Unit.isCoarse(l.unit) ? notCoarse() : l.unit.coarseElements.seq_id_begin.value(l.element)),
    seq_id_end: Element.property(l => !Unit.isCoarse(l.unit) ? notCoarse() : l.unit.coarseElements.seq_id_end.value(l.element)),

    sphere_radius: Element.property(l => !Unit.isSpheres(l.unit) ? notCoarse('spheres') : l.unit.coarseConformation.radius[l.element]),
    sphere_rmsf: Element.property(l => !Unit.isSpheres(l.unit) ? notCoarse('spheres') : l.unit.coarseConformation.rmsf[l.element]),

    gaussian_weight: Element.property(l => !Unit.isGaussians(l.unit) ? notCoarse('gaussians') : l.unit.coarseConformation.weight[l.element]),
    gaussian_covariance_matrix: Element.property(l => !Unit.isGaussians(l.unit) ? notCoarse('gaussians') : l.unit.coarseConformation.covariance_matrix[l.element])
}

function eK(l: Element.Location) {
    switch (l.unit.kind) {
        case Unit.Kind.Atomic:
            return l.unit.model.atomicHierarchy.getEntityKey(l.unit.chainIndex[l.element])
        case Unit.Kind.Spheres:
            return l.unit.model.coarseHierarchy.spheres.entityKey[l.element]
        case Unit.Kind.Gaussians:
            return l.unit.model.coarseHierarchy.gaussians.entityKey[l.element]
    }
}

const entity = {
    key: eK,

    id: Element.property(l => l.unit.model.entities.data.id.value(eK(l))),
    type: Element.property(l => l.unit.model.entities.data.type.value(eK(l))),
    src_method: Element.property(l => l.unit.model.entities.data.src_method.value(eK(l))),
    pdbx_description: Element.property(l => l.unit.model.entities.data.pdbx_description.value(eK(l))),
    formula_weight: Element.property(l => l.unit.model.entities.data.formula_weight.value(eK(l))),
    pdbx_number_of_molecules: Element.property(l => l.unit.model.entities.data.pdbx_number_of_molecules.value(eK(l))),
    details: Element.property(l => l.unit.model.entities.data.details.value(eK(l))),
    pdbx_mutation: Element.property(l => l.unit.model.entities.data.pdbx_mutation.value(eK(l))),
    pdbx_fragment: Element.property(l => l.unit.model.entities.data.pdbx_fragment.value(eK(l))),
    pdbx_ec: Element.property(l => l.unit.model.entities.data.pdbx_ec.value(eK(l)))
}

const unit = {
    operator_name: Element.property(l => l.unit.conformation.operator.name),
    model_num: Element.property(l => l.unit.model.modelNum)
}

const StructureProperties = {
    constant,
    atom,
    residue,
    chain,
    entity,
    unit,
    coarse
}

type StructureProperties = typeof StructureProperties
export default StructureProperties