/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Vec3 } from 'mol-math/linear-algebra';
import { ChunkedArray } from 'mol-data/util';
import { MeshBuilderState } from '../shape/mesh-builder';

const tA = Vec3.zero()
const tB = Vec3.zero()
const tV = Vec3.zero()

const horizontalVector = Vec3.zero()
const verticalVector = Vec3.zero()
const normalOffset = Vec3.zero()
const positionVector = Vec3.zero()
const normalVector = Vec3.zero()
const torsionVector = Vec3.zero()

const p1 = Vec3.zero()
const p2 = Vec3.zero()
const p3 = Vec3.zero()
const p4 = Vec3.zero()

export function addSheet(controlPoints: ArrayLike<number>, normalVectors: ArrayLike<number>, binormalVectors: ArrayLike<number>, linearSegments: number, width: number, height: number, arrowHeight: number, startCap: boolean, endCap: boolean, state: MeshBuilderState) {
    const { vertices, normals, indices } = state

    let vertexCount = vertices.elementCount
    let offsetLength = 0

    if (arrowHeight > 0) {
        Vec3.fromArray(tA, controlPoints, 0)
        Vec3.fromArray(tB, controlPoints, linearSegments * 3)
        offsetLength = arrowHeight / Vec3.magnitude(Vec3.sub(tV, tB, tA))
    }

    for (let i = 0; i <= linearSegments; ++i) {
        const actualHeight = arrowHeight === 0 ? height : arrowHeight * (1 - i / linearSegments);
        const i3 = i * 3

        Vec3.fromArray(verticalVector, normalVectors, i3)
        Vec3.scale(verticalVector, verticalVector, actualHeight);

        Vec3.fromArray(horizontalVector, binormalVectors, i3)
        Vec3.scale(horizontalVector, horizontalVector, width);

        if (arrowHeight > 0) {
            Vec3.fromArray(tA, normalVectors, i3)
            Vec3.fromArray(tB, binormalVectors, i3)
            Vec3.scale(normalOffset, Vec3.cross(normalOffset, tA, tB), offsetLength)
        }

        Vec3.fromArray(positionVector, controlPoints, i3)
        Vec3.fromArray(normalVector, normalVectors, i3)
        Vec3.fromArray(torsionVector, binormalVectors, i3)

        Vec3.add(tA, Vec3.add(tA, positionVector, horizontalVector), verticalVector)
        Vec3.copy(tB, normalVector)
        ChunkedArray.add3(vertices, tA[0], tA[1], tA[2])
        ChunkedArray.add3(normals, tB[0], tB[1], tB[2])

        Vec3.add(tA, Vec3.sub(tA, positionVector, horizontalVector), verticalVector)
        ChunkedArray.add3(vertices, tA[0], tA[1], tA[2])
        ChunkedArray.add3(normals, tB[0], tB[1], tB[2])

        // Vec3.add(tA, Vec3.sub(tA, positionVector, horizontalVector), verticalVector) // reuse tA
        Vec3.add(tB, Vec3.negate(tB, torsionVector), normalOffset)
        ChunkedArray.add3(vertices, tA[0], tA[1], tA[2])
        ChunkedArray.add3(normals, tB[0], tB[1], tB[2])

        Vec3.sub(tA, Vec3.sub(tA, positionVector, horizontalVector), verticalVector)
        ChunkedArray.add3(vertices, tA[0], tA[1], tA[2])
        ChunkedArray.add3(normals, tB[0], tB[1], tB[2])

        // Vec3.sub(tA, Vec3.sub(tA, positionVector, horizontalVector), verticalVector) // reuse tA
        Vec3.negate(tB, normalVector)
        ChunkedArray.add3(vertices, tA[0], tA[1], tA[2])
        ChunkedArray.add3(normals, tB[0], tB[1], tB[2])

        Vec3.sub(tA, Vec3.add(tA, positionVector, horizontalVector), verticalVector)
        ChunkedArray.add3(vertices, tA[0], tA[1], tA[2])
        ChunkedArray.add3(normals, tB[0], tB[1], tB[2])

        // Vec3.sub(tA, Vec3.add(tA, positionVector, horizontalVector), verticalVector) // reuse tA
        Vec3.add(tB, torsionVector, normalOffset)
        ChunkedArray.add3(vertices, tA[0], tA[1], tA[2])
        ChunkedArray.add3(normals, tB[0], tB[1], tB[2])

        Vec3.add(tA, Vec3.add(tA, positionVector, horizontalVector), verticalVector)
        ChunkedArray.add3(vertices, tA[0], tA[1], tA[2])
        ChunkedArray.add3(normals, tB[0], tB[1], tB[2])
    }

    for (let i = 0; i < linearSegments; ++i) {
        for (let j = 0; j < 4; j++) {
            ChunkedArray.add3(
                indices,
                vertexCount + i * 8 + 2 * j,
                vertexCount + (i + 1) * 8 + 2 * j + 1,
                vertexCount + i * 8 + 2 * j + 1
            );
            ChunkedArray.add3(
                indices,
                vertexCount + i * 8 + 2 * j,
                vertexCount + (i + 1) * 8 + 2 * j,
                vertexCount + (i + 1) * 8 + 2 * j + 1
            );
        }
    }

    if (startCap) {
        const offset = 0
        vertexCount = vertices.elementCount

        Vec3.fromArray(verticalVector, normalVectors, offset)
        Vec3.scale(verticalVector, verticalVector, arrowHeight === 0 ? height : arrowHeight);

        Vec3.fromArray(horizontalVector, binormalVectors, offset)
        Vec3.scale(horizontalVector, horizontalVector, width);

        Vec3.fromArray(positionVector, controlPoints, offset)

        Vec3.add(p1, Vec3.add(p1, positionVector, horizontalVector), verticalVector)
        Vec3.sub(p2, Vec3.add(p2, positionVector, horizontalVector), verticalVector)
        Vec3.sub(p3, Vec3.sub(p3, positionVector, horizontalVector), verticalVector)
        Vec3.add(p4, Vec3.sub(p4, positionVector, horizontalVector), verticalVector)

        ChunkedArray.add3(vertices, p1[0], p1[1], p1[2])
        ChunkedArray.add3(vertices, p2[0], p2[1], p2[2])
        ChunkedArray.add3(vertices, p3[0], p3[1], p3[2])
        ChunkedArray.add3(vertices, p4[0], p4[1], p4[2])

        Vec3.cross(normalVector, horizontalVector, verticalVector)

        for (let i = 0; i < 4; ++i) {
            ChunkedArray.add3(normals, normalVector[0], normalVector[1], normalVector[2])
        }
        ChunkedArray.add3(indices, vertexCount + 2, vertexCount + 1, vertexCount);
        ChunkedArray.add3(indices, vertexCount, vertexCount + 3, vertexCount + 2);
    }

    if (endCap && arrowHeight === 0) {
        const offset = linearSegments * 3
        vertexCount = vertices.elementCount

        Vec3.fromArray(verticalVector, normalVectors, offset)
        Vec3.scale(verticalVector, verticalVector, height);

        Vec3.fromArray(horizontalVector, binormalVectors, offset)
        Vec3.scale(horizontalVector, horizontalVector, width);

        Vec3.fromArray(positionVector, controlPoints, offset)

        Vec3.add(p1, Vec3.add(p1, positionVector, horizontalVector), verticalVector)
        Vec3.sub(p2, Vec3.add(p2, positionVector, horizontalVector), verticalVector)
        Vec3.sub(p3, Vec3.sub(p3, positionVector, horizontalVector), verticalVector)
        Vec3.add(p4, Vec3.sub(p4, positionVector, horizontalVector), verticalVector)

        ChunkedArray.add3(vertices, p1[0], p1[1], p1[2])
        ChunkedArray.add3(vertices, p2[0], p2[1], p2[2])
        ChunkedArray.add3(vertices, p3[0], p3[1], p3[2])
        ChunkedArray.add3(vertices, p4[0], p4[1], p4[2])

        Vec3.cross(normalVector, horizontalVector, verticalVector)

        for (let i = 0; i < 4; ++i) {
            ChunkedArray.add3(normals, normalVector[0], normalVector[1], normalVector[2])
        }
        ChunkedArray.add3(indices, vertexCount + 2, vertexCount + 1, vertexCount);
        ChunkedArray.add3(indices, vertexCount, vertexCount + 3, vertexCount + 2);
    }

    return (linearSegments + 1) * 8 + (startCap ? 4 : 0) + (endCap && arrowHeight === 0 ? 4 : 0)
}