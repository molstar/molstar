/**
 * Copyright (c) 2018-2021 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Vec3 } from '../../../../mol-math/linear-algebra';
import { ChunkedArray } from '../../../../mol-data/util';
import { MeshBuilder } from '../mesh-builder';

const tA = Vec3();
const tB = Vec3();
const tV = Vec3();

const horizontalVector = Vec3();
const verticalVector = Vec3();
const verticalRightVector = Vec3();
const verticalLeftVector = Vec3();
const normalOffset = Vec3();
const positionVector = Vec3();
const normalVector = Vec3();
const torsionVector = Vec3();

const p1 = Vec3();
const p2 = Vec3();
const p3 = Vec3();
const p4 = Vec3();

// avoiding namespace lookup improved performance in Chrome (Aug 2020)
const v3fromArray = Vec3.fromArray;
const v3scale = Vec3.scale;
const v3add = Vec3.add;
const v3sub = Vec3.sub;
const v3magnitude = Vec3.magnitude;
const v3negate = Vec3.negate;
const v3copy = Vec3.copy;
const v3cross = Vec3.cross;
const v3set = Vec3.set;
const caAdd3 = ChunkedArray.add3;
const caAdd = ChunkedArray.add;

function addCap(offset: number, state: MeshBuilder.State, controlPoints: ArrayLike<number>, normalVectors: ArrayLike<number>, binormalVectors: ArrayLike<number>, width: number, leftHeight: number, rightHeight: number, flip: boolean) {
    const { vertices, normals, indices } = state;
    const vertexCount = vertices.elementCount;

    v3fromArray(tA, normalVectors, offset);
    v3scale(verticalLeftVector, tA, leftHeight);
    v3scale(verticalRightVector, tA, rightHeight);

    v3fromArray(tB, binormalVectors, offset);
    v3scale(horizontalVector, tB, width);

    v3cross(normalVector, tB, tA);

    v3fromArray(positionVector, controlPoints, offset);

    v3add(p1, v3add(p1, positionVector, horizontalVector), verticalRightVector);
    v3sub(p2, v3add(p2, positionVector, horizontalVector), verticalLeftVector);
    v3sub(p3, v3sub(p3, positionVector, horizontalVector), verticalLeftVector);
    v3add(p4, v3sub(p4, positionVector, horizontalVector), verticalRightVector);

    if (leftHeight < rightHeight) {
        caAdd3(vertices, p4[0], p4[1], p4[2]);
        caAdd3(vertices, p3[0], p3[1], p3[2]);
        caAdd3(vertices, p2[0], p2[1], p2[2]);
        caAdd3(vertices, p1[0], p1[1], p1[2]);
        v3copy(verticalVector, verticalRightVector);
    } else {
        caAdd3(vertices, p1[0], p1[1], p1[2]);
        caAdd3(vertices, p2[0], p2[1], p2[2]);
        caAdd3(vertices, p3[0], p3[1], p3[2]);
        caAdd3(vertices, p4[0], p4[1], p4[2]);
        v3copy(verticalVector, verticalLeftVector);
    }

    if (flip) {
        for (let i = 0; i < 4; ++i) {
            caAdd3(normals, -normalVector[0], -normalVector[1], -normalVector[2]);
        }
        caAdd3(indices, vertexCount, vertexCount + 1, vertexCount + 2);
        caAdd3(indices, vertexCount + 2, vertexCount + 3, vertexCount);
    } else {
        for (let i = 0; i < 4; ++i) {
            caAdd3(normals, normalVector[0], normalVector[1], normalVector[2]);
        }
        caAdd3(indices, vertexCount + 2, vertexCount + 1, vertexCount);
        caAdd3(indices, vertexCount, vertexCount + 3, vertexCount + 2);
    }
}

/** set arrowHeight = 0 for no arrow */
export function addSheet(state: MeshBuilder.State, controlPoints: ArrayLike<number>, normalVectors: ArrayLike<number>, binormalVectors: ArrayLike<number>, linearSegments: number, widthValues: ArrayLike<number>, heightValues: ArrayLike<number>, arrowHeight: number, startCap: boolean, endCap: boolean) {
    const { currentGroup, vertices, normals, indices, groups } = state;

    const vertexCount = vertices.elementCount;
    let offsetLength = 0;

    if (arrowHeight > 0) {
        v3fromArray(tA, controlPoints, 0);
        v3fromArray(tB, controlPoints, linearSegments * 3);
        offsetLength = arrowHeight / v3magnitude(v3sub(tV, tB, tA));
    } else {
        v3set(normalOffset, 0, 0, 0);
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

        v3add(tA, v3add(tA, positionVector, horizontalVector), verticalVector);
        v3add(tB, normalVector, normalOffset);
        caAdd3(vertices, tA[0], tA[1], tA[2]);
        caAdd3(normals, tB[0], tB[1], tB[2]);

        v3add(tA, v3sub(tA, positionVector, horizontalVector), verticalVector);
        caAdd3(vertices, tA[0], tA[1], tA[2]);
        caAdd3(normals, tB[0], tB[1], tB[2]);

        // v3add(tA, v3sub(tA, positionVector, horizontalVector), verticalVector) // reuse tA
        v3negate(tB, torsionVector);
        caAdd3(vertices, tA[0], tA[1], tA[2]);
        caAdd3(normals, tB[0], tB[1], tB[2]);

        v3sub(tA, v3sub(tA, positionVector, horizontalVector), verticalVector);
        caAdd3(vertices, tA[0], tA[1], tA[2]);
        caAdd3(normals, tB[0], tB[1], tB[2]);

        // v3sub(tA, v3sub(tA, positionVector, horizontalVector), verticalVector) // reuse tA
        v3add(tB, v3negate(tB, normalVector), normalOffset);
        caAdd3(vertices, tA[0], tA[1], tA[2]);
        caAdd3(normals, tB[0], tB[1], tB[2]);

        v3sub(tA, v3add(tA, positionVector, horizontalVector), verticalVector);
        caAdd3(vertices, tA[0], tA[1], tA[2]);
        caAdd3(normals, tB[0], tB[1], tB[2]);

        // v3sub(tA, v3add(tA, positionVector, horizontalVector), verticalVector) // reuse tA
        v3copy(tB, torsionVector);
        caAdd3(vertices, tA[0], tA[1], tA[2]);
        caAdd3(normals, tB[0], tB[1], tB[2]);

        v3add(tA, v3add(tA, positionVector, horizontalVector), verticalVector);
        caAdd3(vertices, tA[0], tA[1], tA[2]);
        caAdd3(normals, tB[0], tB[1], tB[2]);
    }

    for (let i = 0; i < linearSegments; ++i) {
        // the triangles are arranged such that opposing triangles of the sheet align
        // which prevents triangle intersection within tight curves
        for (let j = 0; j < 2; j++) {
            caAdd3(
                indices,
                vertexCount + i * 8 + 2 * j, // a
                vertexCount + (i + 1) * 8 + 2 * j + 1, // c
                vertexCount + i * 8 + 2 * j + 1 // b
            );
            caAdd3(
                indices,
                vertexCount + i * 8 + 2 * j, // a
                vertexCount + (i + 1) * 8 + 2 * j, // d
                vertexCount + (i + 1) * 8 + 2 * j + 1 // c
            );
        }
        for (let j = 2; j < 4; j++) {
            caAdd3(
                indices,
                vertexCount + i * 8 + 2 * j, // a
                vertexCount + (i + 1) * 8 + 2 * j, // d
                vertexCount + i * 8 + 2 * j + 1, // b
            );
            caAdd3(
                indices,
                vertexCount + (i + 1) * 8 + 2 * j, // d
                vertexCount + (i + 1) * 8 + 2 * j + 1, // c
                vertexCount + i * 8 + 2 * j + 1, // b
            );
        }
    }

    if (startCap) {
        const width = widthValues[0];
        const height = heightValues[0];
        const h = arrowHeight === 0 ? height : arrowHeight;
        addCap(0, state, controlPoints, normalVectors, binormalVectors, width, h, h, false);
    } else if (arrowHeight > 0) {
        const width = widthValues[0];
        const height = heightValues[0];
        addCap(0, state, controlPoints, normalVectors, binormalVectors, width, arrowHeight, -height, false);
        addCap(0, state, controlPoints, normalVectors, binormalVectors, width, -arrowHeight, height, false);
    }

    if (endCap && arrowHeight === 0) {
        const width = widthValues[linearSegments];
        const height = heightValues[linearSegments];
        addCap(linearSegments * 3, state, controlPoints, normalVectors, binormalVectors, width, height, height, true);
    }

    const addedVertexCount = (linearSegments + 1) * 8 +
        (startCap ? 4 : (arrowHeight > 0 ? 8 : 0)) +
        (endCap && arrowHeight === 0 ? 4 : 0);
    for (let i = 0, il = addedVertexCount; i < il; ++i) caAdd(groups, currentGroup);
}