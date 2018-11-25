/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { AtomicData } from '../atomic';
import { ChemicalComponentMap } from '../chemical-component';
import { AtomicIndex, AtomicDerivedData } from '../atomic/hierarchy';
import { ElementIndex, ResidueIndex } from '../../indexing';
import { MoleculeType } from '../../types';
import { getAtomIdForAtomRole } from 'mol-model/structure/util';

export function getAtomicDerivedData(data: AtomicData, index: AtomicIndex, chemicalComponentMap: ChemicalComponentMap): AtomicDerivedData {
    
    const { label_comp_id, _rowCount: n } = data.residues

    const traceElementIndex = new Uint32Array(n)
    const directionElementIndex = new Uint32Array(n)
    const moleculeType = new Uint8Array(n)

    for (let i = 0; i < n; ++i) {
        const compId = label_comp_id.value(i)
        const chemCompMap = chemicalComponentMap
        const cc = chemCompMap.get(compId)
        const molType = cc ? cc.moleculeType : MoleculeType.unknown
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