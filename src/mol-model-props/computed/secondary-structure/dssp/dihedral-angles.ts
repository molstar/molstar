/**
 * Copyright (c) 2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author Sebastian Bittrich <sebastian.bittrich@rcsb.org>
 */

import { Unit } from '../../../../mol-model/structure';
import { Vec3 } from '../../../../mol-math/linear-algebra';
import { ProteinInfo } from './protein-info';
import { ElementIndex } from '../../../../mol-model/structure/model';
import { radToDeg } from '../../../../mol-math/misc';

export interface DihedralAngles {
    phi: Float32Array
    psi: Float32Array
}

export function calculateUnitDihedralAngles(unit: Unit.Atomic, proteinInfo: ProteinInfo): DihedralAngles {
    const { cIndices, nIndices, residueIndices } = proteinInfo;
    const { position } = unit.conformation;
    const { index } = unit.model.atomicHierarchy;
    const { traceElementIndex } = unit.model.atomicHierarchy.derived.residue;

    const residueCount = residueIndices.length;
    const p = (i: ElementIndex | -1, v: Vec3) => i === -1 ? Vec3.setNaN(v) : position(i, v);

    let cPosPrev = Vec3(), caPosPrev = Vec3(), nPosPrev = Vec3();
    let cPos = Vec3(), caPos = Vec3(), nPos = Vec3();
    let cPosNext = Vec3(), caPosNext = Vec3(), nPosNext = Vec3();

    if (residueCount === 0) return { phi: new Float32Array(0), psi: new Float32Array(0) };

    const phi: Float32Array = new Float32Array(residueCount - 1);
    const psi: Float32Array = new Float32Array(residueCount - 1);

    p(-1, cPosPrev);
    p(-1, caPosPrev);
    p(-1, nPosPrev);

    p(cIndices[0], cPos);
    p(traceElementIndex[residueIndices[0]], caPos);
    p(nIndices[0], nPos);

    p(cIndices[1], cPosNext);
    p(traceElementIndex[residueIndices[1]], caPosNext);
    p(nIndices[1], nPosNext);

    for (let i = 0; i < residueCount - 1; ++i) {
        // ignore C-terminal residue as acceptor
        if (index.findAtomOnResidue(residueIndices[i], 'OXT') !== -1) continue;

        // returns NaN for missing atoms
        phi[i] = radToDeg(Vec3.dihedralAngle(cPosPrev, nPos, caPos, cPos));
        psi[i] = radToDeg(Vec3.dihedralAngle(nPos, caPos, cPos, nPosNext));

        cPosPrev = cPos, caPosPrev = caPos, nPosPrev = nPos;
        cPos = cPosNext, caPos = caPosNext, nPos = nPosNext;

        p(cIndices[i + 1], cPosNext);
        p(traceElementIndex[residueIndices[i + 1]], caPosNext);
        p(nIndices[i + 1], nPosNext);
    }

    return { phi, psi };
}