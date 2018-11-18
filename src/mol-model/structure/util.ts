/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Model, ResidueIndex, ElementIndex } from './model';
import { MoleculeType, AtomRole, MoleculeTypeAtomRoleId } from './model/types';
import { Vec3 } from 'mol-math/linear-algebra';
import { Unit } from './structure';
import Matrix from 'mol-math/linear-algebra/matrix/matrix';

export function getCoarseBegCompId(unit: Unit.Spheres | Unit.Gaussians, element: ElementIndex) {
    const entityKey = unit.coarseElements.entityKey[element]
    const seq = unit.model.sequence.byEntityKey[entityKey]
    const seq_id_begin = unit.coarseElements.seq_id_begin.value(element)
    return seq.compId.value(seq_id_begin - 1) // 1-indexed
}

export function getElementMoleculeType(unit: Unit, element: ElementIndex) {
    let compId = ''
    switch (unit.kind) {
        case Unit.Kind.Atomic:
            compId = unit.model.atomicHierarchy.residues.label_comp_id.value(unit.residueIndex[element])
            break
        case Unit.Kind.Spheres:
        case Unit.Kind.Gaussians:
            compId = getCoarseBegCompId(unit, element)
            break
    }
    const chemCompMap = unit.model.properties.chemicalComponentMap
    const cc = chemCompMap.get(compId)
    return cc ? cc.moleculeType : MoleculeType.unknown
}

export function getAtomicMoleculeType(model: Model, rI: ResidueIndex) {
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

export function getElementIndexForAtomRole(model: Model, rI: ResidueIndex, atomRole: AtomRole) {
    const atomId = getAtomIdForAtomRole(getAtomicMoleculeType(model, rI), atomRole)
    return model.atomicHierarchy.index.findAtomOnResidue(rI, atomId)
}

export function residueLabel(model: Model, rI: number) {
    const { residues, chains, residueAtomSegments, chainAtomSegments } = model.atomicHierarchy
    const { label_comp_id, label_seq_id } = residues
    const { label_asym_id } = chains
    const cI = chainAtomSegments.index[residueAtomSegments.offsets[rI]]
    return `${label_asym_id.value(cI)} ${label_comp_id.value(rI)} ${label_seq_id.value(rI)}`
}

// const centerPos = Vec3.zero()
// const centerMin = Vec3.zero()
// export function getCenterAndRadius(centroid: Vec3, unit: Unit, indices: ArrayLike<number>) {
//     const pos = unit.conformation.position
//     const { elements } = unit
//     Vec3.set(centroid, 0, 0, 0)
//     for (let i = 0, il = indices.length; i < il; ++i) {
//         pos(elements[indices[i]], centerPos)
//         Vec3.add(centroid, centroid, centerPos)
//         Vec3.min(centerMin, centerMin, centerPos)
//     }
//     Vec3.scale(centroid, centroid, 1/indices.length)
//     return Vec3.distance(centerMin, centroid)
// }

const matrixPos = Vec3.zero()
export function getPositionMatrix(unit: Unit, indices: ArrayLike<number>) {
    const pos = unit.conformation.position
    const mat = Matrix.create(3, indices.length)
    const { elements } = unit
    for (let i = 0, il = indices.length; i < il; ++i) {
        pos(elements[indices[i]], matrixPos)
        Vec3.toArray(matrixPos, mat.data, i * 3)
    }
    return mat
}