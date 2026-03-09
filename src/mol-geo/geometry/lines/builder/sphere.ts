/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { LinesBuilder } from '../lines-builder';
import { Vec3 } from '../../../../mol-math/linear-algebra/3d/vec3';
import { Mat4 } from '../../../../mol-math/linear-algebra/3d/mat4';

const _p0 = Vec3();
const _p1 = Vec3();

export interface AddSphereOptions {
    /** Number of segments per circle (default: 32) */
    segments?: number,
    /** Number of circles per dimension, evenly spaced along each axis (default: 1) */
    circlesPerDimension?: number,
}

const DefaultAddSphereOptions: Required<AddSphereOptions> = {
    segments: 32,
    circlesPerDimension: 1,
};

/**
 * Add a wireframe sphere to a LinesBuilder as orthogonal circles in XY, XZ, YZ planes.
 */
export function addSphere(builder: LinesBuilder, radius: number, transform: Mat4, group: number, options?: AddSphereOptions) {
    const segments = options?.segments ?? DefaultAddSphereOptions.segments;
    const circlesPerDim = options?.circlesPerDimension ?? DefaultAddSphereOptions.circlesPerDimension;

    // For each dimension, draw `circlesPerDim` circles spaced evenly along the axis.
    // circlesPerDim=1 → one circle at the center (offset=0)
    // circlesPerDim=3 → circles at -r/2, 0, +r/2

    for (let dim = 0; dim < 3; ++dim) {
        for (let ci = 0; ci < circlesPerDim; ++ci) {
            // offset along the perpendicular axis: evenly spaced in [-radius, radius]
            const offset = circlesPerDim === 1
                ? 0
                : -radius + (2 * radius * (ci + 1)) / (circlesPerDim + 1);

            // Choose a smaller radius for offset circles (cross-section of the sphere)
            const r = Math.sqrt(Math.max(0, radius * radius - offset * offset));
            if (r < 1e-8) continue;

            addCircle(builder, r, offset, dim, transform, group, segments);
        }
    }
}

function addCircle(builder: LinesBuilder, radius: number, offset: number, perpAxis: number, transform: Mat4, group: number, segments: number) {
    // perpAxis: 0=X (circle in YZ), 1=Y (circle in XZ), 2=Z (circle in XY)
    for (let i = 0; i < segments; ++i) {
        const a0 = (i / segments) * Math.PI * 2;
        const a1 = ((i + 1) / segments) * Math.PI * 2;

        const cos0 = Math.cos(a0) * radius;
        const sin0 = Math.sin(a0) * radius;
        const cos1 = Math.cos(a1) * radius;
        const sin1 = Math.sin(a1) * radius;

        if (perpAxis === 2) {
            // XY circle at z=offset
            Vec3.set(_p0, cos0, sin0, offset);
            Vec3.set(_p1, cos1, sin1, offset);
        } else if (perpAxis === 1) {
            // XZ circle at y=offset
            Vec3.set(_p0, cos0, offset, sin0);
            Vec3.set(_p1, cos1, offset, sin1);
        } else {
            // YZ circle at x=offset
            Vec3.set(_p0, offset, cos0, sin0);
            Vec3.set(_p1, offset, cos1, sin1);
        }

        Vec3.transformMat4(_p0, _p0, transform);
        Vec3.transformMat4(_p1, _p1, transform);
        builder.addVec(_p0, _p1, group);
    }
}
