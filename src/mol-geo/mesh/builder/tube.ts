/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Vec3 } from 'mol-math/linear-algebra';
import { ChunkedArray } from 'mol-data/util';
import { MeshBuilder } from '../mesh-builder';

const normalVector = Vec3.zero()
const binormalVector = Vec3.zero()
const controlPoint = Vec3.zero()
const tempPos = Vec3.zero()
const a = Vec3.zero()
const b = Vec3.zero()
const u = Vec3.zero()
const v = Vec3.zero()

export function addTube(builder: MeshBuilder, controlPoints: ArrayLike<number>, normalVectors: ArrayLike<number>, binormalVectors: ArrayLike<number>, linearSegments: number, radialSegments: number, width: number, height: number, waveFactor: number, startCap: boolean, endCap: boolean) {
    const { currentGroup, vertices, normals, indices, groups } = builder.state

    let vertexCount = vertices.elementCount
    const di = 1 / linearSegments

    for (let i = 0; i <= linearSegments; ++i) {
        const i3 = i * 3
        Vec3.fromArray(u, normalVectors, i3)
        Vec3.fromArray(v, binormalVectors, i3)

        const tt = di * i - 0.5;
        const ff = 1 + (waveFactor - 1) * (Math.cos(2 * Math.PI * tt) + 1);
        const w = ff * width, h = ff * height;

        for (let j = 0; j < radialSegments; ++j) {
            const t = 2 * Math.PI * j / radialSegments;

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
        }
    }

    for (let i = 0; i < linearSegments; ++i) {
        for (let j = 0; j < radialSegments; ++j) {
            ChunkedArray.add3(
                indices,
                vertexCount + i * radialSegments + (j + 1) % radialSegments,
                vertexCount + (i + 1) * radialSegments + (j + 1) % radialSegments,
                vertexCount + i * radialSegments + j
            );
            ChunkedArray.add3(
                indices,
                vertexCount + (i + 1) * radialSegments + (j + 1) % radialSegments,
                vertexCount + (i + 1) * radialSegments + j,
                vertexCount + i * radialSegments + j
            );
        }
    }

    if (startCap) {
        const offset = 0
        const centerVertex = vertices.elementCount
        Vec3.fromArray(u, normalVectors, offset)
        Vec3.fromArray(v, binormalVectors, offset)
        Vec3.fromArray(controlPoint, controlPoints, offset)
        Vec3.cross(normalVector, u, v)

        ChunkedArray.add3(vertices, controlPoint[0], controlPoint[1], controlPoint[2]);
        ChunkedArray.add3(normals, normalVector[0], normalVector[1], normalVector[2]);

        vertexCount = vertices.elementCount
        for (let i = 0; i < radialSegments; ++i) {
            const t = 2 * Math.PI * i / radialSegments;

            Vec3.copy(a, u)
            Vec3.copy(b, v)
            Vec3.add(
                tempPos,
                Vec3.scale(a, a, height * Math.cos(t)),
                Vec3.scale(b, b, width * Math.sin(t))
            )
            Vec3.add(tempPos, controlPoint, tempPos)

            ChunkedArray.add3(vertices, tempPos[0], tempPos[1], tempPos[2]);
            ChunkedArray.add3(normals, normalVector[0], normalVector[1], normalVector[2]);

            ChunkedArray.add3(
                indices,
                centerVertex,
                vertexCount + i,
                vertexCount + (i + 1) % radialSegments
            );
        }
    }

    if (endCap) {
        const offset = linearSegments * 3
        const centerVertex = vertices.elementCount
        Vec3.fromArray(u, normalVectors, offset)
        Vec3.fromArray(v, binormalVectors, offset)
        Vec3.fromArray(controlPoint, controlPoints, offset)
        Vec3.cross(normalVector, u, v)

        ChunkedArray.add3(vertices, controlPoint[0], controlPoint[1], controlPoint[2]);
        ChunkedArray.add3(normals, normalVector[0], normalVector[1], normalVector[2]);

        vertexCount = vertices.elementCount
        for (let i = 0; i < radialSegments; ++i) {
            const t = 2 * Math.PI * i / radialSegments

            Vec3.copy(a, u)
            Vec3.copy(b, v)
            Vec3.add(
                tempPos,
                Vec3.scale(a, a, height * Math.cos(t)),
                Vec3.scale(b, b, width * Math.sin(t))
            )
            Vec3.add(tempPos, controlPoint, tempPos)

            ChunkedArray.add3(vertices, tempPos[0], tempPos[1], tempPos[2]);
            ChunkedArray.add3(normals, normalVector[0], normalVector[1], normalVector[2]);

            ChunkedArray.add3(
                indices,
                vertexCount + i,
                vertexCount + (i + 1) % radialSegments,
                centerVertex
            );
        }
    }

    const addedVertexCount = (linearSegments + 1) * radialSegments + (startCap ? radialSegments + 1 : 0) + (endCap ? radialSegments + 1 : 0)
    for (let i = 0, il = addedVertexCount; i < il; ++i) ChunkedArray.add(groups, currentGroup)
}