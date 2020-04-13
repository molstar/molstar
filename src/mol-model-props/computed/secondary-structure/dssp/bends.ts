/**
 * Copyright (c) 2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author Sebastian Bittrich <sebastian.bittrich@rcsb.org>
 */

import { Vec3 } from '../../../../mol-math/linear-algebra';
import { radToDeg } from '../../../../mol-math/misc';
import { DSSPContext, DSSPType } from './common';
import { ElementIndex } from '../../../../mol-model/structure';

/**
 * Bend(i) =: [angle ((CW - Ca(i - 2)),(C"(i + 2) - C"(i))) > 70"]
 *
 * Type: S
 */
export function assignBends(ctx: DSSPContext) {
    const { unit, flags, proteinInfo } = ctx;
    const { position } = unit.conformation;
    const { traceElementIndex } = unit.model.atomicHierarchy.derived.residue;

    const { residueIndices, nIndices } = proteinInfo;
    const residueCount = residueIndices.length;

    // const position = (i: number, v: Vec3) => Vec3.set(v, x[i], y[i], z[i])
    const p = (i: ElementIndex | -1, v: Vec3) => i === -1 ? Vec3.setNaN(v) : position(i, v);

    const caPosPrev2 = Vec3();
    const caPos = Vec3();
    const caPosNext2 = Vec3();

    const cPos = Vec3();
    const nPosNext = Vec3();

    const caMinus2 = Vec3();
    const caPlus2 = Vec3();

    f1: for (let i = 2; i < residueCount - 2; i++) {
        // check for peptide bond
        for (let k = 0; k < 4; k++) {
            let index = i + k - 2;
            p(traceElementIndex[index], cPos);
            p(nIndices[index + 1], nPosNext);
            if (Vec3.squaredDistance(cPos, nPosNext) > 6.25 /* max squared peptide bond distance allowed */) {
                continue f1;
            }
        }

        const oRIprev2 = residueIndices[i - 2];
        const oRI = residueIndices[i];
        const oRInext2 = residueIndices[i + 2];

        const caAtomPrev2 = traceElementIndex[oRIprev2];
        const caAtom = traceElementIndex[oRI];
        const caAtomNext2 = traceElementIndex[oRInext2];

        p(caAtomPrev2, caPosPrev2);
        p(caAtom, caPos);
        p(caAtomNext2, caPosNext2);

        Vec3.sub(caMinus2, caPosPrev2, caPos);
        Vec3.sub(caPlus2, caPos, caPosNext2);

        const angle = radToDeg(Vec3.angle(caMinus2, caPlus2));
        if (angle && angle > 70.00) {
            flags[i] |= DSSPType.Flag.S;
        }
    }
}