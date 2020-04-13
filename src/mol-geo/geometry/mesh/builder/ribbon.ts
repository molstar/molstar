/**
 * Copyright (c) 2018-2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Vec3 } from '../../../../mol-math/linear-algebra';
import { ChunkedArray } from '../../../../mol-data/util';
import { MeshBuilder } from '../mesh-builder';

const tA = Vec3.zero();
const tB = Vec3.zero();
const tV = Vec3.zero();

const horizontalVector = Vec3.zero();
const verticalVector = Vec3.zero();
const normalOffset = Vec3.zero();
const positionVector = Vec3.zero();
const normalVector = Vec3.zero();
const torsionVector = Vec3.zero();

/** set arrowHeight = 0 for no arrow */
export function addRibbon(state: MeshBuilder.State, controlPoints: ArrayLike<number>, normalVectors: ArrayLike<number>, binormalVectors: ArrayLike<number>, linearSegments: number, widthValues: ArrayLike<number>, heightValues: ArrayLike<number>, arrowHeight: number) {
    const { currentGroup, vertices, normals, indices, groups } = state;

    let vertexCount = vertices.elementCount;
    let offsetLength = 0;

    if (arrowHeight > 0) {
        Vec3.fromArray(tA, controlPoints, 0);
        Vec3.fromArray(tB, controlPoints, linearSegments * 3);
        offsetLength = arrowHeight / Vec3.magnitude(Vec3.sub(tV, tB, tA));
    }

    for (let i = 0; i <= linearSegments; ++i) {
        const width = widthValues[i];
        const height = heightValues[i];

        const actualHeight = arrowHeight === 0 ? height : arrowHeight * (1 - i / linearSegments);
        const i3 = i * 3;

        Vec3.fromArray(verticalVector, normalVectors, i3);
        Vec3.scale(verticalVector, verticalVector, actualHeight);

        Vec3.fromArray(horizontalVector, binormalVectors, i3);
        Vec3.scale(horizontalVector, horizontalVector, width);

        if (arrowHeight > 0) {
            Vec3.fromArray(tA, normalVectors, i3);
            Vec3.fromArray(tB, binormalVectors, i3);
            Vec3.scale(normalOffset, Vec3.cross(normalOffset, tA, tB), offsetLength);
        }

        Vec3.fromArray(positionVector, controlPoints, i3);
        Vec3.fromArray(normalVector, normalVectors, i3);
        Vec3.fromArray(torsionVector, binormalVectors, i3);

        Vec3.add(tA, positionVector, verticalVector);
        Vec3.negate(tB, torsionVector);
        ChunkedArray.add3(vertices, tA[0], tA[1], tA[2]);
        ChunkedArray.add3(normals, tB[0], tB[1], tB[2]);

        Vec3.sub(tA, positionVector, verticalVector);
        ChunkedArray.add3(vertices, tA[0], tA[1], tA[2]);
        ChunkedArray.add3(normals, tB[0], tB[1], tB[2]);

        Vec3.add(tA, positionVector, verticalVector);
        Vec3.copy(tB, torsionVector);
        ChunkedArray.add3(vertices, tA[0], tA[1], tA[2]);
        ChunkedArray.add3(normals, tB[0], tB[1], tB[2]);

        Vec3.sub(tA, positionVector, verticalVector);
        ChunkedArray.add3(vertices, tA[0], tA[1], tA[2]);
        ChunkedArray.add3(normals, tB[0], tB[1], tB[2]);
    }

    for (let i = 0; i < linearSegments; ++i) {
        ChunkedArray.add3(
            indices,
            vertexCount + i * 4,
            vertexCount + (i + 1) * 4 + 1,
            vertexCount + i * 4 + 1
        );
        ChunkedArray.add3(
            indices,
            vertexCount + i * 4,
            vertexCount + (i + 1) * 4,
            vertexCount + (i + 1) * 4 + 1
        );

        ChunkedArray.add3(
            indices,
            vertexCount + i * 4 + 2 + 1,
            vertexCount + (i + 1) * 4 + 2 + 1,
            vertexCount + i * 4 + 2
        );
        ChunkedArray.add3(
            indices,
            vertexCount + i * 4 + 2,
            vertexCount + (i + 1) * 4 + 2 + 1,
            vertexCount + (i + 1) * 4 + 2
        );
    }

    const addedVertexCount = (linearSegments + 1) * 4;
    for (let i = 0, il = addedVertexCount; i < il; ++i) ChunkedArray.add(groups, currentGroup);
}