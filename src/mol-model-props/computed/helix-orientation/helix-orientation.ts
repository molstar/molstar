/**
 * Copyright (c) 2021 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { ElementIndex } from '../../../mol-model/structure';
import { Segmentation } from '../../../mol-data/int/segmentation';
import { SortedRanges } from '../../../mol-data/int/sorted-ranges';
import { OrderedSet } from '../../../mol-data/int';
import { Model } from '../../../mol-model/structure/model';
import { Vec3 } from '../../../mol-math/linear-algebra';

export interface HelixOrientation {
    centers: ArrayLike<number>
}

/** Usees same definition as GROMACS' helixorient */
export function calcHelixOrientation(model: Model): HelixOrientation {
    const { x, y, z } = model.atomicConformation;
    const { polymerType, traceElementIndex } = model.atomicHierarchy.derived.residue;
    const n = polymerType.length;

    const elements = OrderedSet.ofBounds(0, model.atomicConformation.atomId.rowCount) as OrderedSet<ElementIndex>;
    const polymerIt = SortedRanges.transientSegments(model.atomicRanges.polymerRanges, elements);
    const residueIt = Segmentation.transientSegments(model.atomicHierarchy.residueAtomSegments, elements);

    const centers = new Float32Array(n * 3);
    const axes = new Float32Array(n * 3);

    let i = 0;
    let j = -1;
    let s = -1;

    const a1 = Vec3();
    const a2 = Vec3();
    const a3 = Vec3();
    const a4 = Vec3();

    const r12 = Vec3();
    const r23 = Vec3();
    const r34 = Vec3();

    const v1 = Vec3();
    const v2 = Vec3();
    const vt = Vec3();

    const diff13 = Vec3();
    const diff24 = Vec3();

    const axis = Vec3();
    const prevAxis = Vec3();

    while (polymerIt.hasNext) {
        const ps = polymerIt.move();
        residueIt.setSegment(ps);
        i = -1;
        s = -1;
        while (residueIt.hasNext) {
            i += 1;
            const { index } = residueIt.move();
            if (i === 0) s = index;

            j = (index - 2);
            const j3 = j * 3;

            Vec3.copy(a1, a2);
            Vec3.copy(a2, a3);
            Vec3.copy(a3, a4);

            const eI = traceElementIndex[index];
            Vec3.set(a4, x[eI], y[eI], z[eI]);

            if (i < 3) continue;

            Vec3.sub(r12, a2, a1);
            Vec3.sub(r23, a3, a2);
            Vec3.sub(r34, a4, a3);

            Vec3.sub(diff13, r12, r23);
            Vec3.sub(diff24, r23, r34);

            Vec3.cross(axis, diff13, diff24);
            Vec3.normalize(axis, axis);
            Vec3.toArray(axis, axes, j3);

            const tmp = Math.cos(Vec3.angle(diff13, diff24));

            const diff13Length = Vec3.magnitude(diff13);
            const diff24Length = Vec3.magnitude(diff24);

            const r = (
                Math.sqrt(diff24Length * diff13Length) /
                // clamp, to avoid numerical instabilities for when
                // angle between diff13 and diff24 is close to 0
                Math.max(2.0, 2.0 * (1.0 - tmp))
            );

            Vec3.scale(v1, diff13, r / diff13Length);
            Vec3.sub(v1, a2, v1);
            Vec3.toArray(v1, centers, j3);

            Vec3.scale(v2, diff24, r / diff24Length);
            Vec3.sub(v2, a3, v2);
            Vec3.toArray(v2, centers, j3 + 3);

            Vec3.copy(prevAxis, axis);
        }

        // calc axis as dir of second and third center pos
        // project first trace atom onto axis to get first center pos
        const s3 = s * 3;
        Vec3.fromArray(v1, centers, s3 + 3);
        Vec3.fromArray(v2, centers, s3 + 6);
        Vec3.normalize(axis, Vec3.sub(axis, v1, v2));
        const sI = traceElementIndex[s];
        Vec3.set(a1, x[sI], y[sI], z[sI]);
        Vec3.copy(vt, a1);
        Vec3.projectPointOnVector(vt, vt, axis, v1);
        Vec3.toArray(vt, centers, s3);

        // calc axis as dir of n-1 and n-2 center pos
        // project last traceAtom onto axis to get last center pos
        const e = j + 2;
        const e3 = e * 3;
        Vec3.fromArray(v1, centers, e3 - 3);
        Vec3.fromArray(v2, centers, e3 - 6);
        Vec3.normalize(axis, Vec3.sub(axis, v1, v2));
        const eI = traceElementIndex[e];
        Vec3.set(a1, x[eI], y[eI], z[eI]);
        Vec3.copy(vt, a1);
        Vec3.projectPointOnVector(vt, vt, axis, v1);
        Vec3.toArray(vt, centers, e3);
    }

    return {
        centers
    };
}
