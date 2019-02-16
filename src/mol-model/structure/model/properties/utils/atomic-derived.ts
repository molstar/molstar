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
import { getAtomIdForAtomRole } from 'mol-model/structure/util';

export function getAtomicDerivedData(data: AtomicData, index: AtomicIndex, chemicalComponentMap: ChemicalComponentMap): AtomicDerivedData {
    const { label_comp_id, _rowCount: n } = data.residues

    const traceElementIndex = new Uint32Array(n)
    const directionElementIndex = new Uint32Array(n)
    const moleculeType = new Uint8Array(n)

    const moleculeTypeMap = new Map<string, MoleculeType>()

    for (let i = 0; i < n; ++i) {
        const compId = label_comp_id.value(i)
        const chemCompMap = chemicalComponentMap
        let molType: MoleculeType
        if (moleculeTypeMap.has(compId)){
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
        traceElementIndex[i] = index.findAtomOnResidue(i as ResidueIndex, traceAtomId)

        const directionAtomId = getAtomIdForAtomRole(molType, 'direction')
        directionElementIndex[i] = index.findAtomOnResidue(i as ResidueIndex, directionAtomId)
    }

    return {
        residue: {
            traceElementIndex: traceElementIndex as unknown as ArrayLike<ElementIndex>,
            directionElementIndex: directionElementIndex as unknown as ArrayLike<ElementIndex>,
            moleculeType: moleculeType as unknown as ArrayLike<MoleculeType>,
        }
    }
}