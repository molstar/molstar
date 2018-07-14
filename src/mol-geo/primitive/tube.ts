/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Vec3 } from 'mol-math/linear-algebra';
import { ChunkedArray } from 'mol-data/util';

const normalVector = Vec3.zero()
const binormalVector = Vec3.zero()
const tempPos = Vec3.zero()
const a = Vec3.zero()
const b = Vec3.zero()
const u = Vec3.zero()
const v = Vec3.zero()

export function addTube(controlPoints: Helpers.NumberArray, normalVectors: Helpers.NumberArray, binormalVectors: Helpers.NumberArray, linearSegments: number, radialSegments: number, width: number, height: number, waveFactor: number, vertices: ChunkedArray<number, 3>, normals: ChunkedArray<number, 3>, indices: ChunkedArray<number, 3>, ids: ChunkedArray<number, 1>, currentId: number) {
    const vertexCount = vertices.elementCount
    const di = 1 / linearSegments

    for (let i = 0; i <= linearSegments; ++i) {
        const i3 = i * 3
        Vec3.fromArray(u, normalVectors, i3)
        Vec3.fromArray(v, binormalVectors, i3)

        const tt = di * i - 0.5;
        const ff = 1 + (waveFactor - 1) * (Math.cos(2 * Math.PI * tt) + 1);
        const w = ff * width, h = ff * height;

        for (let j = 0; j < radialSegments; ++j) {
            let t = 2 * Math.PI * j / radialSegments;

            Vec3.copy(a, u)
            Vec3.copy(b, v)
            Vec3.add(
                normalVector,
                Vec3.scale(a, a, w * Math.cos(t)),
                Vec3.scale(b, b, h * Math.sin(t))
            )

            Vec3.copy(a, u)
            Vec3.copy(b, v)
            Vec3.add(
                binormalVector,
                Vec3.scale(a, a, h * Math.cos(t)),
                Vec3.scale(b, b, w * Math.sin(t))
            )
            Vec3.normalize(normalVector, normalVector)

            Vec3.fromArray(tempPos, controlPoints, i3)
            Vec3.add(tempPos, tempPos, binormalVector)

            ChunkedArray.add3(vertices, tempPos[0], tempPos[1], tempPos[2]);
            ChunkedArray.add3(normals, normalVector[0], normalVector[1], normalVector[2]);
            ChunkedArray.add(ids, currentId);
        }
    }

    for (let i = 0; i < linearSegments; ++i) {
        for (let j = 0; j < radialSegments; ++j) {
            ChunkedArray.add3(
                indices,
                (vertexCount + i * radialSegments + (j + 1) % radialSegments),
                (vertexCount + (i + 1) * radialSegments + (j + 1) % radialSegments),
                (vertexCount + i * radialSegments + j)
            );
            ChunkedArray.add3(
                indices,
                (vertexCount + (i + 1) * radialSegments + (j + 1) % radialSegments),
                (vertexCount + (i + 1) * radialSegments + j),
                (vertexCount + i * radialSegments + j)
            );
        }
    }
}