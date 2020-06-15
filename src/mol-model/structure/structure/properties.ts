/**
 * Copyright (c) 2017-2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import StructureElement from './element';
import Unit from './unit';
import { VdwRadius } from '../model/properties/atomic';
import { SecondaryStructureType } from '../model/types';
import { SecondaryStructureProvider } from '../../../mol-model-props/computed/secondary-structure';
import { SymmetryOperator } from '../../../mol-math/geometry';

function p<T>(p: StructureElement.Property<T>) { return p; }

const constant = {
    true: p(l => true),
    false: p(l => false),
    zero: p(l => 0)
};

function notAtomic(): never {
    throw 'Property only available for atomic models.';
}

function notCoarse(kind?: string): never {
    if (!!kind) throw `Property only available for coarse models (${kind}).`;
    throw `Property only available for coarse models.`;
}

// TODO: remove the type checks?

const atom = {
    key: p(l => l.element),

    // Conformation
    x: p(l => l.unit.conformation.x(l.element)),
    y: p(l => l.unit.conformation.y(l.element)),
    z: p(l => l.unit.conformation.z(l.element)),
    id: p(l => !Unit.isAtomic(l.unit) ? notAtomic() : l.unit.model.atomicConformation.atomId.value(l.element)),
    occupancy: p(l => !Unit.isAtomic(l.unit) ?  notAtomic() : l.unit.model.atomicConformation.occupancy.value(l.element)),
    B_iso_or_equiv: p(l => !Unit.isAtomic(l.unit) ?  notAtomic() : l.unit.model.atomicConformation.B_iso_or_equiv.value(l.element)),
    sourceIndex: p(l => Unit.isAtomic(l.unit)
        ? l.unit.model.atomicHierarchy.atoms.sourceIndex.value(l.element)
        // TODO: when implemented, this should map to the source index.
        : l.element),

    // Hierarchy
    type_symbol: p(l => !Unit.isAtomic(l.unit) ? notAtomic() : l.unit.model.atomicHierarchy.atoms.type_symbol.value(l.element)),
    label_atom_id: p(l => !Unit.isAtomic(l.unit) ? notAtomic() : l.unit.model.atomicHierarchy.atoms.label_atom_id.value(l.element)),
    auth_atom_id: p(l => !Unit.isAtomic(l.unit) ? notAtomic() : l.unit.model.atomicHierarchy.atoms.auth_atom_id.value(l.element)),
    label_alt_id: p(l => !Unit.isAtomic(l.unit) ? notAtomic() : l.unit.model.atomicHierarchy.atoms.label_alt_id.value(l.element)),
    label_comp_id: p(compId),
    auth_comp_id: p(l => !Unit.isAtomic(l.unit) ? notAtomic() : l.unit.model.atomicHierarchy.atoms.auth_comp_id.value(l.unit.residueIndex[l.element])),
    pdbx_formal_charge: p(l => !Unit.isAtomic(l.unit) ? notAtomic() : l.unit.model.atomicHierarchy.atoms.pdbx_formal_charge.value(l.element)),

    // Derived
    vdw_radius: p(l => !Unit.isAtomic(l.unit) ? notAtomic() : VdwRadius(l.unit.model.atomicHierarchy.atoms.type_symbol.value(l.element))),
};

function compId(l: StructureElement.Location) {
    if (!Unit.isAtomic(l.unit)) notAtomic();
    return l.unit.model.atomicHierarchy.atoms.label_comp_id.value(l.element);
}

function seqId(l: StructureElement.Location) {
    return !Unit.isAtomic(l.unit) ? notAtomic() : l.unit.model.atomicHierarchy.residues.label_seq_id.value(l.unit.residueIndex[l.element]);
}

function hasMicroheterogeneity(l: StructureElement.Location) {
    if (!Unit.isAtomic(l.unit)) notAtomic();
    const entitySeq = l.unit.model.sequence.byEntityKey[eK(l)];
    return entitySeq && entitySeq.sequence.microHet.has(seqId(l));
}

function microheterogeneityCompIds(l: StructureElement.Location) {
    if (!Unit.isAtomic(l.unit)) notAtomic();
    const entitySeq = l.unit.model.sequence.byEntityKey[eK(l)];
    if (entitySeq) {
        return entitySeq.sequence.microHet.get(seqId(l)) || [compId(l)];
    } else {
        return [compId(l)];
    }
}

const residue = {
    key: p(l => !Unit.isAtomic(l.unit) ? notAtomic() : l.unit.residueIndex[l.element]),

    group_PDB: p(l => !Unit.isAtomic(l.unit) ? notAtomic() : l.unit.model.atomicHierarchy.residues.group_PDB.value(l.unit.residueIndex[l.element])),
    label_seq_id: p(seqId),
    auth_seq_id: p(l => !Unit.isAtomic(l.unit) ? notAtomic() : l.unit.model.atomicHierarchy.residues.auth_seq_id.value(l.unit.residueIndex[l.element])),
    pdbx_PDB_ins_code: p(l => !Unit.isAtomic(l.unit) ? notAtomic() : l.unit.model.atomicHierarchy.residues.pdbx_PDB_ins_code.value(l.unit.residueIndex[l.element])),

    // Properties
    isNonStandard: p(l => !Unit.isAtomic(l.unit) ? notAtomic() : microheterogeneityCompIds(l).some(c => l.unit.model.properties.chemicalComponentMap.get(c)!.mon_nstd_flag[0] !== 'y')),
    hasMicroheterogeneity: p(hasMicroheterogeneity),
    microheterogeneityCompIds: p(microheterogeneityCompIds),
    secondary_structure_type: p(l => {
        if (!Unit.isAtomic(l.unit)) notAtomic();
        const secStruc = SecondaryStructureProvider.get(l.structure).value?.get(l.unit.id);
        return secStruc?.type[l.unit.residueIndex[l.element]] ?? SecondaryStructureType.Flag.NA;
    }),
    secondary_structure_key: p(l => {
        if (!Unit.isAtomic(l.unit)) notAtomic();
        const secStruc = SecondaryStructureProvider.get(l.structure).value?.get(l.unit.id);
        return secStruc?.key[l.unit.residueIndex[l.element]] ?? -1;
    }),
    chem_comp_type: p(l => !Unit.isAtomic(l.unit) ? notAtomic() : l.unit.model.properties.chemicalComponentMap.get(compId(l))!.type),
};

const chain = {
    key: p(l => !Unit.isAtomic(l.unit) ? notAtomic() : l.unit.chainIndex[l.element]),

    label_asym_id: p(l => !Unit.isAtomic(l.unit) ? notAtomic() : l.unit.model.atomicHierarchy.chains.label_asym_id.value(l.unit.chainIndex[l.element])),
    auth_asym_id: p(l => !Unit.isAtomic(l.unit) ? notAtomic() : l.unit.model.atomicHierarchy.chains.auth_asym_id.value(l.unit.chainIndex[l.element])),
    label_entity_id: p(l => !Unit.isAtomic(l.unit) ? notAtomic() : l.unit.model.atomicHierarchy.chains.label_entity_id.value(l.unit.chainIndex[l.element]))
};

const coarse = {
    key: atom.key,
    entityKey: p(l => !Unit.isCoarse(l.unit) ? notCoarse() : l.unit.coarseElements.entityKey[l.element]),

    x: atom.x,
    y: atom.y,
    z: atom.z,

    asym_id: p(l => !Unit.isCoarse(l.unit) ? notCoarse() : l.unit.coarseElements.asym_id.value(l.element)),
    seq_id_begin: p(l => !Unit.isCoarse(l.unit) ? notCoarse() : l.unit.coarseElements.seq_id_begin.value(l.element)),
    seq_id_end: p(l => !Unit.isCoarse(l.unit) ? notCoarse() : l.unit.coarseElements.seq_id_end.value(l.element)),

    sphere_radius: p(l => !Unit.isSpheres(l.unit) ? notCoarse('spheres') : l.unit.coarseConformation.radius[l.element]),
    sphere_rmsf: p(l => !Unit.isSpheres(l.unit) ? notCoarse('spheres') : l.unit.coarseConformation.rmsf[l.element]),

    gaussian_weight: p(l => !Unit.isGaussians(l.unit) ? notCoarse('gaussians') : l.unit.coarseConformation.weight[l.element]),
    gaussian_covariance_matrix: p(l => !Unit.isGaussians(l.unit) ? notCoarse('gaussians') : l.unit.coarseConformation.covariance_matrix[l.element])
};

function eK(l: StructureElement.Location) {
    switch (l.unit.kind) {
        case Unit.Kind.Atomic:
            return l.unit.model.atomicHierarchy.index.getEntityFromChain(l.unit.chainIndex[l.element]);
        case Unit.Kind.Spheres:
            return l.unit.model.coarseHierarchy.spheres.entityKey[l.element];
        case Unit.Kind.Gaussians:
            return l.unit.model.coarseHierarchy.gaussians.entityKey[l.element];
    }
}

const entity = {
    key: p(eK),

    id: p(l => l.unit.model.entities.data.id.value(eK(l))),
    type: p(l => l.unit.model.entities.data.type.value(eK(l))),
    src_method: p(l => l.unit.model.entities.data.src_method.value(eK(l))),
    pdbx_description: p(l => l.unit.model.entities.data.pdbx_description.value(eK(l))),
    formula_weight: p(l => l.unit.model.entities.data.formula_weight.value(eK(l))),
    pdbx_number_of_molecules: p(l => l.unit.model.entities.data.pdbx_number_of_molecules.value(eK(l))),
    details: p(l => l.unit.model.entities.data.details.value(eK(l))),
    pdbx_mutation: p(l => l.unit.model.entities.data.pdbx_mutation.value(eK(l))),
    pdbx_fragment: p(l => l.unit.model.entities.data.pdbx_fragment.value(eK(l))),
    pdbx_ec: p(l => l.unit.model.entities.data.pdbx_ec.value(eK(l))),

    subtype: p(l => l.unit.model.entities.subtype.value(eK(l))),
    prd_id: p(l => l.unit.model.entities.prd_id.value(eK(l))),
};

const _emptyList: any[] = [];
const unit = {
    id: p(l => l.unit.id),
    chainGroupId: p(l => l.unit.chainGroupId),
    multiChain: p(l => Unit.Traits.is(l.unit.traits, Unit.Trait.MultiChain)),
    object_primitive: p(l => l.unit.objectPrimitive),
    operator_name: p(l => l.unit.conformation.operator.name),
    model_index: p(l => l.unit.model.modelNum),
    model_label: p(l => l.unit.model.label),
    model_entry_id: p(l => l.unit.model.entryId),
    hkl: p(l => l.unit.conformation.operator.hkl),
    spgrOp: p(l => l.unit.conformation.operator.spgrOp),

    model_num: p(l => l.unit.model.modelNum),
    pdbx_struct_assembly_id: p(l => l.unit.conformation.operator.assembly?.id || SymmetryOperator.DefaultName),
    pdbx_struct_oper_list_ids: p(l => l.unit.conformation.operator.assembly?.operList || _emptyList),
    struct_ncs_oper_id: p(l => l.unit.conformation.operator.ncsId),
};

const StructureProperties = {
    constant,
    atom,
    residue,
    chain,
    entity,
    unit,
    coarse
};

type StructureProperties = typeof StructureProperties
export default StructureProperties;