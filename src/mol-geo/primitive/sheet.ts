/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Vec3 } from 'mol-math/linear-algebra';
import { ChunkedArray } from 'mol-data/util';

const tA = Vec3.zero()
const tB = Vec3.zero()
const tV = Vec3.zero()

const horizontalVector = Vec3.zero()
const verticalVector = Vec3.zero()
const normalOffset = Vec3.zero()
const positionVector = Vec3.zero()
const normalVector = Vec3.zero()
const torsionVector = Vec3.zero()

export function addSheet(controlPoints: Helpers.NumberArray, normalVectors: Helpers.NumberArray, binormalVectors: Helpers.NumberArray, linearSegments: number, width: number, height: number, arrowWidth: number, vertices: ChunkedArray<number, 3>, normals: ChunkedArray<number, 3>, indices: ChunkedArray<number, 3>, ids: ChunkedArray<number, 1>, currentId: number) {

    const vertexCount = vertices.elementCount
    let offsetLength = 0

    if (arrowWidth > 0) {
        Vec3.fromArray(tA, controlPoints, 0)
        Vec3.fromArray(tB, controlPoints, linearSegments * 3)
        offsetLength = arrowWidth / Vec3.magnitude(Vec3.sub(tV, tB, tA))
    }

    for (let i = 0; i <= linearSegments; ++i) {
        const actualWidth = arrowWidth === 0 ? width : arrowWidth * (1 - i / linearSegments);

        const i3 = i * 3

        Vec3.fromArray(verticalVector, normalVectors, i3)
        Vec3.scale(verticalVector, verticalVector, actualWidth);

        Vec3.fromArray(horizontalVector, binormalVectors, i3)
        Vec3.scale(horizontalVector, horizontalVector, height);

        if (arrowWidth > 0) {
            Vec3.fromArray(tA, normalVectors, i3)
            Vec3.fromArray(tB, binormalVectors, i3)
            Vec3.scale(normalOffset, Vec3.cross(normalOffset, tA, tB), offsetLength)
        }

        Vec3.fromArray(positionVector, controlPoints, i3)
        Vec3.fromArray(normalVector, normalVectors, i3)
        Vec3.fromArray(torsionVector, binormalVectors, i3)

        Vec3.add(tA, Vec3.add(tA, Vec3.copy(tA, positionVector), horizontalVector), verticalVector)
        Vec3.copy(tB, normalVector)
        ChunkedArray.add3(vertices, tA[0], tA[1], tA[2])
        ChunkedArray.add3(normals, tB[0], tB[1], tB[2])
        ChunkedArray.add(ids, currentId);

        Vec3.add(tA, Vec3.sub(tA, Vec3.copy(tA, positionVector), horizontalVector), verticalVector)
        ChunkedArray.add3(vertices, tA[0], tA[1], tA[2])
        ChunkedArray.add3(normals, tB[0], tB[1], tB[2])
        ChunkedArray.add(ids, currentId);

        Vec3.add(tA, Vec3.sub(tA, Vec3.copy(tA, positionVector), horizontalVector), verticalVector)
        Vec3.add(tB, Vec3.scale(tB, Vec3.copy(tB, torsionVector), -1), normalOffset)
        ChunkedArray.add3(vertices, tA[0], tA[1], tA[2])
        ChunkedArray.add3(normals, tB[0], tB[1], tB[2])
        ChunkedArray.add(ids, currentId);

        Vec3.sub(tA, Vec3.sub(tA, Vec3.copy(tA, positionVector), horizontalVector), verticalVector)
        ChunkedArray.add3(vertices, tA[0], tA[1], tA[2])
        ChunkedArray.add3(normals, tB[0], tB[1], tB[2])
        ChunkedArray.add(ids, currentId);

        Vec3.sub(tA, Vec3.sub(tA, Vec3.copy(tA, positionVector), horizontalVector), verticalVector)
        Vec3.scale(tB, Vec3.copy(tB, normalVector), -1)
        ChunkedArray.add3(vertices, tA[0], tA[1], tA[2])
        ChunkedArray.add3(normals, tB[0], tB[1], tB[2])
        ChunkedArray.add(ids, currentId);

        Vec3.sub(tA, Vec3.add(tA, Vec3.copy(tA, positionVector), horizontalVector), verticalVector)
        ChunkedArray.add3(vertices, tA[0], tA[1], tA[2])
        ChunkedArray.add3(normals, tB[0], tB[1], tB[2])
        ChunkedArray.add(ids, currentId);

        Vec3.sub(tA, Vec3.add(tA, Vec3.copy(tA, positionVector), horizontalVector), verticalVector)
        Vec3.add(tB, Vec3.copy(tB, torsionVector), normalOffset)
        ChunkedArray.add3(vertices, tA[0], tA[1], tA[2])
        ChunkedArray.add3(normals, tB[0], tB[1], tB[2])
        ChunkedArray.add(ids, currentId);

        Vec3.add(tA, Vec3.add(tA, Vec3.copy(tA, positionVector), horizontalVector), verticalVector)
        ChunkedArray.add3(vertices, tA[0], tA[1], tA[2])
        ChunkedArray.add3(normals, tB[0], tB[1], tB[2])
        ChunkedArray.add(ids, currentId);
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
}