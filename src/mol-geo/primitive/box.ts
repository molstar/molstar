/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

// adapted from three.js, MIT License Copyright 2010-2018 three.js authors

import { Vec3 } from 'mol-math/linear-algebra'

export const DefaultBoxProps = {
    width: 1,
    height: 1,
    depth: 1,
    widthSegments: 1,
    heightSegments: 1,
    depthSegments: 1
}
export type BoxProps = typeof DefaultBoxProps

export default function Box(props?: Partial<BoxProps>) {
    const { width, height, depth, widthSegments, heightSegments, depthSegments } = { ...DefaultBoxProps, ...props }

    // buffers
    const indices: number[] = [];
    const vertices: number[] = [];
    const normals: number[] = [];

    // helper variables
    let numberOfVertices = 0;

    // build each side of the box geometry
    buildPlane(2, 1, 0, -1, -1, depth, height, width, depthSegments, heightSegments); // px
    buildPlane(2, 1, 0, 1, -1, depth, height, -width, depthSegments, heightSegments); // nx
    buildPlane(0, 2, 1, 1, 1, width, depth, height, widthSegments, depthSegments); // py
    buildPlane(0, 2, 1, 1, -1, width, depth, -height, widthSegments, depthSegments); // ny
    buildPlane(0, 1, 2, 1, -1, width, height, depth, widthSegments, heightSegments); // pz
    buildPlane(0, 1, 2, -1, -1, width, height, -depth, widthSegments, heightSegments); // nz

    return {
        vertices: new Float32Array(vertices),
        normals: new Float32Array(normals),
        indices: new Uint32Array(indices)
    }

    function buildPlane(u: number, v: number, w: number, udir: number, vdir: number, width: number, height: number, depth: number, gridX: number, gridY: number) {

        const segmentWidth = width / gridX;
        const segmentHeight = height / gridY;

        const widthHalf = width / 2;
        const heightHalf = height / 2;
        const depthHalf = depth / 2;

        const gridX1 = gridX + 1;
        const gridY1 = gridY + 1;

        let vertexCounter = 0;

        const vector = Vec3.zero();

        // generate vertices and normals
        for (let iy = 0; iy < gridY1; ++iy) {
            const y = iy * segmentHeight - heightHalf;
            for (let ix = 0; ix < gridX1; ++ix) {
                const x = ix * segmentWidth - widthHalf;

                // set values to correct vector component
                vector[ u ] = x * udir;
                vector[ v ] = y * vdir;
                vector[ w ] = depthHalf;

                // now apply vector to vertex buffer
                vertices.push(...vector);

                // set values to correct vector component
                vector[ u ] = 0;
                vector[ v ] = 0;
                vector[ w ] = depth > 0 ? 1 : -1;

                // now apply vector to normal buffer
                normals.push(...vector);

                vertexCounter += 1;
            }
        }

        // indices
        // 1. you need three indices to draw a single face
        // 2. a single segment consists of two faces
        // 3. so we need to generate six (2*3) indices per segment
        for (let iy = 0; iy < gridY; ++iy) {
            for (let ix = 0; ix < gridX; ++ix) {
                const a = numberOfVertices + ix + gridX1 * iy;
                const b = numberOfVertices + ix + gridX1 * (iy + 1);
                const c = numberOfVertices + (ix + 1) + gridX1 * (iy + 1);
                const d = numberOfVertices + (ix + 1) + gridX1 * iy;

                // faces
                indices.push(a, b, d);
                indices.push(b, c, d);
            }
        }

        numberOfVertices += vertexCounter;
    }
}