/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { LinesBuilder } from '../lines-builder';
import { Vec3 } from '../../../../mol-math/linear-algebra/3d/vec3';
import { Mat4 } from '../../../../mol-math/linear-algebra/3d/mat4';

const _start = Vec3();
const _end = Vec3();

// 8 corners of the unit cube [0,1]^3: indexed by 3 bits (x=bit0, y=bit1, z=bit2)
const _corners: Vec3[] = [];
for (let i = 0; i < 8; ++i) _corners.push(Vec3());

// 12 box edges: [startCorner, endCorner, axis]
// axis: 0=x, 1=y, 2=z — the axis along which the edge runs
const BoxEdges: [number, number, number][] = [
    // X-axis edges (y,z vary)
    [0, 1, 0], [2, 3, 0], [4, 5, 0], [6, 7, 0],
    // Y-axis edges (x,z vary)
    [0, 2, 1], [1, 3, 1], [4, 6, 1], [5, 7, 1],
    // Z-axis edges (x,y vary)
    [0, 4, 2], [1, 5, 2], [2, 6, 2], [3, 7, 2],
];

/**
 * Add a wireframe box to a LinesBuilder.
 */
export function addBox(builder: LinesBuilder, transform: Mat4, group: number) {
    // Compute 8 corners in world space from unit cube [0,1]^3
    // Corner index bits: bit0=x(0/1), bit1=y(0/1), bit2=z(0/1)
    for (let ci = 0; ci < 8; ++ci) {
        Vec3.set(_corners[ci],
            (ci & 1) ? 1 : 0,
            (ci & 2) ? 1 : 0,
            (ci & 4) ? 1 : 0,
        );
        Vec3.transformMat4(_corners[ci], _corners[ci], transform);
    }

    for (const [si, ei] of BoxEdges) {
        Vec3.copy(_start, _corners[si]);
        Vec3.copy(_end, _corners[ei]);

        // Draw edge line
        builder.addVec(_start, _end, group);
    }
}
