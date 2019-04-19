/**
 * Copyright (c) 2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author Sebastian Bittrich <sebastian.bittrich@rcsb.org>
 */

import { AtomicHierarchy, AtomicConformation } from 'mol-model/structure/model/properties/atomic';
import { BackboneAtomIndices } from './common';
import { ResidueIndex } from 'mol-model/structure';
import { SortedArray } from 'mol-data/int';
import { Vec3 } from 'mol-math/linear-algebra';

export function calculateDihedralAngles(hierarchy: AtomicHierarchy, conformation: AtomicConformation, proteinResidues: SortedArray<ResidueIndex>, backboneIndices: BackboneAtomIndices): { phi: Float32Array, psi: Float32Array } {
    const { cIndices, nIndices } = backboneIndices
    const { index } = hierarchy
    const { x, y, z } = conformation
    const { traceElementIndex } = hierarchy.derived.residue

    const residueCount = proteinResidues.length
    const position = (i: number, v: Vec3) => i === -1 ? Vec3.setNaN(v) : Vec3.set(v, x[i], y[i], z[i])

    let cPosPrev = Vec3(), caPosPrev = Vec3(), nPosPrev = Vec3()
    let cPos = Vec3(), caPos = Vec3(), nPos = Vec3()
    let cPosNext = Vec3(), caPosNext = Vec3(), nPosNext = Vec3()

    if (residueCount === 0) return { phi: new Float32Array(0), psi: new Float32Array(0) }

    const phi: Float32Array = new Float32Array(residueCount - 1)
    const psi: Float32Array = new Float32Array(residueCount - 1)

    position(-1, cPosPrev)
    position(-1, caPosPrev)
    position(-1, nPosPrev)

    position(cIndices[0], cPos)
    position(traceElementIndex[proteinResidues[0]], caPos)
    position(nIndices[0], nPos)

    position(cIndices[1], cPosNext)
    position(traceElementIndex[proteinResidues[1]], caPosNext)
    position(nIndices[1], nPosNext)

    for (let i = 0; i < residueCount - 1; ++i) {
        // ignore C-terminal residue as acceptor
        if (index.findAtomOnResidue(proteinResidues[i], 'OXT') !== -1) continue

        // returns NaN for missing atoms
        phi[i] = Vec3.dihedralAngle(cPosPrev, nPos, caPos, cPos)
        psi[i] = Vec3.dihedralAngle(nPos, caPos, cPos, nPosNext)

        cPosPrev = cPos, caPosPrev = caPos, nPosPrev = nPos
        cPos = cPosNext, caPos = caPosNext, nPos = nPosNext

        position(cIndices[i + 1], cPosNext)
        position(traceElementIndex[proteinResidues[i + 1]], caPosNext)
        position(nIndices[i + 1], nPosNext)
    }

    return { phi, psi };
}