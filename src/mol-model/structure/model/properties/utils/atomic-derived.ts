/**
 * Copyright (c) 2018-2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { AtomicData } from '../atomic';
import { ChemicalComponentMap } from '../chemical-component';
import { AtomicIndex, AtomicDerivedData } from '../atomic/hierarchy';
import { ElementIndex, ResidueIndex } from '../../indexing';
import { MoleculeType, getMoleculeType, getComponentType } from '../../types';
import { getAtomIdForAtomRole } from '../../../../../mol-model/structure/util';

export function getAtomicDerivedData(data: AtomicData, index: AtomicIndex, chemicalComponentMap: ChemicalComponentMap): AtomicDerivedData {
    const { label_comp_id, _rowCount: n } = data.residues

    const traceElementIndex = new Int32Array(n)
    const directionFromElementIndex = new Int32Array(n)
    const directionToElementIndex = new Int32Array(n)
    const moleculeType = new Uint8Array(n)

    const moleculeTypeMap = new Map<string, MoleculeType>()

    for (let i = 0; i < n; ++i) {
        const compId = label_comp_id.value(i)
        const chemCompMap = chemicalComponentMap
        let molType: MoleculeType
        if (moleculeTypeMap.has(compId)) {
            molType = moleculeTypeMap.get(compId)!
        } else if (chemCompMap.has(compId)) {
            molType = getMoleculeType(chemCompMap.get(compId)!.type, compId)
            moleculeTypeMap.set(compId, molType)
        } else {
            molType = getMoleculeType(getComponentType(compId), compId)
            // TODO if unknown molecule type, use atom names to guess molecule type
            moleculeTypeMap.set(compId, molType)
        }
        moleculeType[i] = molType

        const traceAtomId = getAtomIdForAtomRole(molType, 'trace')
        let traceIndex = index.findAtomsOnResidue(i as ResidueIndex, traceAtomId)
        if (traceIndex === -1) {
            const coarseAtomId = getAtomIdForAtomRole(molType, 'coarseBackbone')
            traceIndex = index.findAtomsOnResidue(i as ResidueIndex, coarseAtomId)
        }
        traceElementIndex[i] = traceIndex

        const directionFromAtomId = getAtomIdForAtomRole(molType, 'directionFrom')
        directionFromElementIndex[i] = index.findAtomsOnResidue(i as ResidueIndex, directionFromAtomId)

        const directionToAtomId = getAtomIdForAtomRole(molType, 'directionTo')
        directionToElementIndex[i] = index.findAtomsOnResidue(i as ResidueIndex, directionToAtomId)
    }

    return {
        residue: {
            traceElementIndex: traceElementIndex as unknown as ArrayLike<ElementIndex | -1>,
            directionFromElementIndex: directionFromElementIndex as unknown as ArrayLike<ElementIndex | -1>,
            directionToElementIndex: directionToElementIndex as unknown as ArrayLike<ElementIndex | -1>,
            moleculeType: moleculeType as unknown as ArrayLike<MoleculeType>,
        }
    }
}