/**
 * Copyright (c) 2018-2025 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Primitive } from './primitive';
import { Cage } from './cage';

const plane: Primitive = {
    vertices: new Float32Array([
        -0.5, 0.5, 0,
        0.5, 0.5, 0,
        -0.5, -0.5, 0,
        0.5, -0.5, 0
    ]),
    normals: new Float32Array([
        0, 0, 1,
        0, 0, 1,
        0, 0, 1,
        0, 0, 1
    ]),
    indices: new Uint32Array([
        0, 2, 1,
        1, 2, 3
    ])
};

const planeCage: Cage = {
    vertices: plane.vertices,
    edges: new Uint32Array([0, 1, 2, 3, 3, 1, 2, 0])
};

export function Plane(): Primitive {
    return plane;
}

export function PlaneCage(): Cage {
    return planeCage;
}

//

export function SegmentedPlane(widthSegments: number, heightSegments: number): Primitive {
    const widthSegments1 = widthSegments + 1;
    const heightSegments1 = heightSegments + 1;

    const segmentWidth = 1 / widthSegments;
    const segmentHeight = 1 / heightSegments;

    const vertices = new Float32Array(widthSegments1 * heightSegments1 * 3);
    const normals = new Float32Array(widthSegments1 * heightSegments1 * 3);
    const indices = new Uint32Array(widthSegments * heightSegments * 6);

    let i = 0;
    for (let iy = 0; iy < heightSegments1; ++iy) {
        const y = iy * segmentHeight - 0.5;
        for (let ix = 0; ix < widthSegments1; ++ix) {
            const x = ix * segmentWidth - 0.5;
            vertices[i] = x;
            vertices[i + 1] = -y;
            vertices[i + 2] = 0;
            normals[i] = 0;
            normals[i + 1] = 0;
            normals[i + 2] = 1;
            i += 3;
        }
    }

    let j = 0;
    for (let iy = 0; iy < heightSegments; ++iy) {
        for (let ix = 0; ix < widthSegments; ++ix) {
            const a = ix + widthSegments1 * iy;
            const b = ix + widthSegments1 * (iy + 1);
            const c = (ix + 1) + widthSegments1 * (iy + 1);
            const d = (ix + 1) + widthSegments1 * iy;

            indices[j] = a;
            indices[j + 1] = b;
            indices[j + 2] = d;
            indices[j + 3] = b;
            indices[j + 4] = c;
            indices[j + 5] = d;
            j += 6;
        }
    }

    return { vertices, normals, indices };
}
