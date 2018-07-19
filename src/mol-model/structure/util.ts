/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Model, ResidueIndex, ElementIndex } from './model';
import { MoleculeType } from './model/types';

export function getMoleculeType(model: Model, rI: ResidueIndex) {
    const compId = model.atomicHierarchy.residues.label_comp_id.value(rI)
    const chemCompMap = model.properties.chemicalComponentMap
    const cc = chemCompMap.get(compId)
    return cc ? cc.moleculeType : MoleculeType.unknown
}

export function getElementIndexForAtomId(model: Model, rI: ResidueIndex, atomId: string): ElementIndex {
    const { offsets } = model.atomicHierarchy.residueAtomSegments
    const { label_atom_id } = model.atomicHierarchy.atoms
    for (let j = offsets[rI], _j = offsets[rI + 1]; j < _j; j++) {
        if (label_atom_id.value(j) === atomId) return j as ElementIndex
    }
    return offsets[rI] as ElementIndex
}

export function residueLabel(model: Model, rI: number) {
    const { residues, chains, residueAtomSegments, chainAtomSegments } = model.atomicHierarchy
    const { label_comp_id, label_seq_id } = residues
    const { label_asym_id } = chains
    const cI = chainAtomSegments.index[residueAtomSegments.offsets[rI]]
    return `${label_asym_id.value(cI)} ${label_comp_id.value(rI)} ${label_seq_id.value(rI)}`
}