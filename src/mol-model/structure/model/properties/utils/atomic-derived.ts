/**
 * Copyright (c) 2018-2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { AtomicData } from '../atomic';
import { AtomicIndex, AtomicDerivedData, AtomicSegments } from '../atomic/hierarchy';
import { ElementIndex, ResidueIndex } from '../../indexing';
import { MoleculeType, getMoleculeType, getComponentType, PolymerType, getPolymerType } from '../../types';
import { getAtomIdForAtomRole } from '../../../../../mol-model/structure/util';
import { ChemicalComponentMap } from '../common';
import { isProductionMode } from '../../../../../mol-util/debug';

export function getAtomicDerivedData(data: AtomicData, segments: AtomicSegments, index: AtomicIndex, chemicalComponentMap: ChemicalComponentMap): AtomicDerivedData {
    const { label_comp_id } = data.atoms;
    const { _rowCount: n } = data.residues;
    const { offsets } = segments.residueAtomSegments;

    const traceElementIndex = new Int32Array(n);
    const directionFromElementIndex = new Int32Array(n);
    const directionToElementIndex = new Int32Array(n);
    const moleculeType = new Uint8Array(n);
    const polymerType = new Uint8Array(n);

    const moleculeTypeMap = new Map<string, MoleculeType>();
    const polymerTypeMap = new Map<string, PolymerType>();

    for (let i = 0 as ResidueIndex; i < n; ++i) {
        const compId = label_comp_id.value(offsets[i]);
        const chemCompMap = chemicalComponentMap;

        let molType: MoleculeType;
        let polyType: PolymerType;
        if (moleculeTypeMap.has(compId)) {
            molType = moleculeTypeMap.get(compId)!;
            polyType = polymerTypeMap.get(compId)!;
        } else {
            let type: string;
            if (chemCompMap.has(compId)) {
                type = chemCompMap.get(compId)!.type;
            } else {
                if (!isProductionMode) console.info('chemComp not found', compId);
                type = getComponentType(compId);
            }
            molType = getMoleculeType(type, compId);
            // TODO if unknown molecule type, use atom names to guess molecule type
            polyType = getPolymerType(type, molType);
            moleculeTypeMap.set(compId, molType);
            polymerTypeMap.set(compId, polyType);
        }
        moleculeType[i] = molType;
        polymerType[i] = polyType;

        const traceAtomId = getAtomIdForAtomRole(polyType, 'trace');
        let traceIndex = index.findAtomsOnResidue(i, traceAtomId);
        if (traceIndex === -1) {
            const coarseAtomId = getAtomIdForAtomRole(polyType, 'coarseBackbone');
            traceIndex = index.findAtomsOnResidue(i, coarseAtomId);
        }
        traceElementIndex[i] = traceIndex;

        const directionFromAtomId = getAtomIdForAtomRole(polyType, 'directionFrom');
        directionFromElementIndex[i] = index.findAtomsOnResidue(i, directionFromAtomId);

        const directionToAtomId = getAtomIdForAtomRole(polyType, 'directionTo');
        directionToElementIndex[i] = index.findAtomsOnResidue(i, directionToAtomId);
    }

    return {
        residue: {
            traceElementIndex: traceElementIndex as unknown as ArrayLike<ElementIndex | -1>,
            directionFromElementIndex: directionFromElementIndex as unknown as ArrayLike<ElementIndex | -1>,
            directionToElementIndex: directionToElementIndex as unknown as ArrayLike<ElementIndex | -1>,
            moleculeType: moleculeType as unknown as ArrayLike<MoleculeType>,
            polymerType: polymerType as unknown as ArrayLike<PolymerType>,
        }
    };
}