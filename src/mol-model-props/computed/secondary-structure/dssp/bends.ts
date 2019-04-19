/**
 * Copyright (c) 2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author Sebastian Bittrich <sebastian.bittrich@rcsb.org>
 */

import { Vec3 } from 'mol-math/linear-algebra';
import { radToDeg } from 'mol-math/misc';
import { DSSPContext, DSSPType } from './common';

/**
 * Bend(i) =: [angle ((CW - Ca(i - 2)),(C"(i + 2) - C"(i))) > 70"]
 *
 * Type: S
 */
export function assignBends(ctx: DSSPContext) {
    const flags = ctx.flags
    const { x, y, z } = ctx.conformation
    const { traceElementIndex } = ctx.hierarchy.derived.residue

    const proteinResidues = ctx.proteinResidues
    const residueCount = proteinResidues.length

    const position = (i: number, v: Vec3) => Vec3.set(v, x[i], y[i], z[i])

    const caPosPrev2 = Vec3()
    const caPos = Vec3()
    const caPosNext2 = Vec3()

    const nIndices = ctx.backboneIndices.nIndices
    const cPos = Vec3()
    const nPosNext = Vec3()

    const caMinus2 = Vec3()
    const caPlus2 = Vec3()

    f1: for (let i = 2; i < residueCount - 2; i++) {
        // check for peptide bond
        for (let k = 0; k < 4; k++) {
            let index = i + k - 2
            position(traceElementIndex[index], cPos)
            position(nIndices[index + 1], nPosNext)
            if (Vec3.squaredDistance(cPos, nPosNext) > 6.25 /* max squared peptide bond distance allowed */) {
                continue f1
            }
        }

        const oRIprev2 = proteinResidues[i - 2]
        const oRI = proteinResidues[i]
        const oRInext2 = proteinResidues[i + 2]

        const caAtomPrev2 = traceElementIndex[oRIprev2]
        const caAtom = traceElementIndex[oRI]
        const caAtomNext2 = traceElementIndex[oRInext2]

        position(caAtomPrev2, caPosPrev2)
        position(caAtom, caPos)
        position(caAtomNext2, caPosNext2)

        Vec3.sub(caMinus2, caPosPrev2, caPos)
        Vec3.sub(caPlus2, caPos, caPosNext2)

        const angle = radToDeg(Vec3.angle(caMinus2, caPlus2))
        if (angle && angle > 70.00) {
            flags[i] |= DSSPType.Flag.S
        }
    }
}