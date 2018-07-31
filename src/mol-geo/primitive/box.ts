/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

// adapted from three.js, MIT License Copyright 2010-2018 three.js authors

import { Vec3 } from 'mol-math/linear-algebra'
import { Primitive } from './primitive';

export const DefaultBoxProps = {
    width: 1,
    height: 1,
    depth: 1
}
export type BoxProps = Partial<typeof DefaultBoxProps>

const tmpVector = Vec3.zero();

export function Box(props?: BoxProps): Primitive {
    const { width, height, depth } = { ...DefaultBoxProps, ...props }

    // buffers
    const vertices = new Float32Array(72);
    const normals = new Float32Array(72);
    const indices = new Uint32Array(36);

    // helper variables
    let vertexCount = 0;

    // build each side of the box geometry
    buildPlane(2, 1, 0, -1, -1, depth, height, width); // px
    buildPlane(2, 1, 0, 1, -1, depth, height, -width); // nx
    buildPlane(0, 2, 1, 1, 1, width, depth, height); // py
    buildPlane(0, 2, 1, 1, -1, width, depth, -height); // ny
    buildPlane(0, 1, 2, 1, -1, width, height, depth); // pz
    buildPlane(0, 1, 2, -1, -1, width, height, -depth); // nz

    return { vertices, normals, indices }

    function buildPlane(u: number, v: number, w: number, udir: number, vdir: number, width: number, height: number, depth: number) {
        // generate vertices and normals
        for (let iy = 0; iy < 2; ++iy) {
            const y = iy * height - height / 2;
            for (let ix = 0; ix < 2; ++ix) {
                const x = ix * width - width / 2;

                // set values to correct vector component and add to vertex buffer
                tmpVector[u] = x * udir;
                tmpVector[v] = y * vdir;
                tmpVector[w] = depth / 2;
                Vec3.toArray(tmpVector, vertices, vertexCount * 3);

                // set values to correct vector component and add to normal buffer
                tmpVector[u] = 0;
                tmpVector[v] = 0;
                tmpVector[w] = depth > 0 ? 1 : -1;
                Vec3.toArray(tmpVector, normals, vertexCount * 3);

                ++vertexCount;
            }
        }

        // faces
        const vc = vertexCount - 4
        const iidx = (vc / 2) * 3
        indices[iidx] = vc
        indices[iidx + 1] = vc + 2
        indices[iidx + 2] = vc + 1
        indices[iidx + 3] = vc + 2
        indices[iidx + 4] = vc + 3
        indices[iidx + 5] = vc + 1
    }
}