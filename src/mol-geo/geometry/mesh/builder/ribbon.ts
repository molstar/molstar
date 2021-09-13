/**
 * Copyright (c) 2018-2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Vec3 } from '../../../../mol-math/linear-algebra';
import { ChunkedArray } from '../../../../mol-data/util';
import { MeshBuilder } from '../mesh-builder';

// avoiding namespace lookup improved performance in Chrome (Aug 2020)
const v3fromArray = Vec3.fromArray;
const v3magnitude = Vec3.magnitude;
const v3sub = Vec3.sub;
const v3add = Vec3.add;
const v3scale = Vec3.scale;
const v3negate = Vec3.negate;
const v3copy = Vec3.copy;
const v3cross = Vec3.cross;
const caAdd3 = ChunkedArray.add3;
const caAdd = ChunkedArray.add;

const tA = Vec3();
const tB = Vec3();
const tV = Vec3();

const horizontalVector = Vec3();
const verticalVector = Vec3();
const normalOffset = Vec3();
const positionVector = Vec3();
const normalVector = Vec3();
const torsionVector = Vec3();

/** set arrowHeight = 0 for no arrow */
export function addRibbon(state: MeshBuilder.State, controlPoints: ArrayLike<number>, normalVectors: ArrayLike<number>, binormalVectors: ArrayLike<number>, linearSegments: number, widthValues: ArrayLike<number>, heightValues: ArrayLike<number>, arrowHeight: number) {
    const { currentGroup, vertices, normals, indices, groups } = state;

    const vertexCount = vertices.elementCount;
    let offsetLength = 0;

    if (arrowHeight > 0) {
        v3fromArray(tA, controlPoints, 0);
        v3fromArray(tB, controlPoints, linearSegments * 3);
        offsetLength = arrowHeight / v3magnitude(v3sub(tV, tB, tA));
    }

    for (let i = 0; i <= linearSegments; ++i) {
        const width = widthValues[i];
        const height = heightValues[i];

        const actualHeight = arrowHeight === 0 ? height : arrowHeight * (1 - i / linearSegments);
        const i3 = i * 3;

        v3fromArray(verticalVector, normalVectors, i3);
        v3scale(verticalVector, verticalVector, actualHeight);

        v3fromArray(horizontalVector, binormalVectors, i3);
        v3scale(horizontalVector, horizontalVector, width);

        if (arrowHeight > 0) {
            v3fromArray(tA, normalVectors, i3);
            v3fromArray(tB, binormalVectors, i3);
            v3scale(normalOffset, v3cross(normalOffset, tA, tB), offsetLength);
        }

        v3fromArray(positionVector, controlPoints, i3);
        v3fromArray(normalVector, normalVectors, i3);
        v3fromArray(torsionVector, binormalVectors, i3);

        v3add(tA, positionVector, verticalVector);
        v3negate(tB, torsionVector);
        caAdd3(vertices, tA[0], tA[1], tA[2]);
        caAdd3(normals, tB[0], tB[1], tB[2]);

        v3sub(tA, positionVector, verticalVector);
        caAdd3(vertices, tA[0], tA[1], tA[2]);
        caAdd3(normals, tB[0], tB[1], tB[2]);

        v3add(tA, positionVector, verticalVector);
        v3copy(tB, torsionVector);
        caAdd3(vertices, tA[0], tA[1], tA[2]);
        caAdd3(normals, tB[0], tB[1], tB[2]);

        v3sub(tA, positionVector, verticalVector);
        caAdd3(vertices, tA[0], tA[1], tA[2]);
        caAdd3(normals, tB[0], tB[1], tB[2]);
    }

    for (let i = 0; i < linearSegments; ++i) {
        caAdd3(
            indices,
            vertexCount + i * 4,
            vertexCount + (i + 1) * 4 + 1,
            vertexCount + i * 4 + 1
        );
        caAdd3(
            indices,
            vertexCount + i * 4,
            vertexCount + (i + 1) * 4,
            vertexCount + (i + 1) * 4 + 1
        );

        caAdd3(
            indices,
            vertexCount + i * 4 + 2 + 1,
            vertexCount + (i + 1) * 4 + 2 + 1,
            vertexCount + i * 4 + 2
        );
        caAdd3(
            indices,
            vertexCount + i * 4 + 2,
            vertexCount + (i + 1) * 4 + 2 + 1,
            vertexCount + (i + 1) * 4 + 2
        );
    }

    const addedVertexCount = (linearSegments + 1) * 4;
    for (let i = 0, il = addedVertexCount; i < il; ++i) caAdd(groups, currentGroup);
}