/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Model, ResidueIndex, ElementIndex } from './model';
import { MoleculeType, AtomRole, MoleculeTypeAtomRoleId } from './model/types';
import { Vec3 } from 'mol-math/linear-algebra';
import { Unit } from './structure';

export function getMoleculeType(model: Model, rI: ResidueIndex) {
    const compId = model.atomicHierarchy.residues.label_comp_id.value(rI)
    const chemCompMap = model.properties.chemicalComponentMap
    const cc = chemCompMap.get(compId)
    return cc ? cc.moleculeType : MoleculeType.unknown
}

export function getAtomIdForAtomRole(moleculeType: MoleculeType, atomRole: AtomRole) {
    const m = MoleculeTypeAtomRoleId[moleculeType]
    if (m !== undefined) {
        const a = m[atomRole]
        if (a !== undefined) return a
    }
    return ''
}

export function getElementIndexForAtomId(model: Model, rI: ResidueIndex, atomId: string): ElementIndex {
    const { offsets } = model.atomicHierarchy.residueAtomSegments
    const { label_atom_id } = model.atomicHierarchy.atoms
    for (let j = offsets[rI], _j = offsets[rI + 1]; j < _j; j++) {
        if (label_atom_id.value(j) === atomId) return j as ElementIndex
    }
    return offsets[rI] as ElementIndex
}

export function getElementIndexForAtomRole(model: Model, rI: ResidueIndex, atomRole: AtomRole) {
    const atomId = getAtomIdForAtomRole(getMoleculeType(model, rI), atomRole)
    return getElementIndexForAtomId(model, rI, atomId)
}

export function residueLabel(model: Model, rI: number) {
    const { residues, chains, residueAtomSegments, chainAtomSegments } = model.atomicHierarchy
    const { label_comp_id, label_seq_id } = residues
    const { label_asym_id } = chains
    const cI = chainAtomSegments.index[residueAtomSegments.offsets[rI]]
    return `${label_asym_id.value(cI)} ${label_comp_id.value(rI)} ${label_seq_id.value(rI)}`
}

const centerPos = Vec3.zero()
const centerMin = Vec3.zero()
export function getCenterAndRadius(centroid: Vec3, unit: Unit, indices: ArrayLike<number>) {
    const pos = unit.conformation.position
    const { elements } = unit
    Vec3.set(centroid, 0, 0, 0)
    for (let i = 0, il = indices.length; i < il; ++i) {
        pos(elements[indices[i]], centerPos)
        Vec3.add(centroid, centroid, centerPos)
        Vec3.min(centerMin, centerMin, centerPos)
    }
    Vec3.scale(centroid, centroid, 1/indices.length)
    return Vec3.distance(centerMin, centroid)
}