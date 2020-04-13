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
const verticalRightVector = Vec3.zero();
const verticalLeftVector = Vec3.zero();
const normalOffset = Vec3.zero();
const positionVector = Vec3.zero();
const normalVector = Vec3.zero();
const torsionVector = Vec3.zero();

const p1 = Vec3.zero();
const p2 = Vec3.zero();
const p3 = Vec3.zero();
const p4 = Vec3.zero();

function addCap(offset: number, state: MeshBuilder.State, controlPoints: ArrayLike<number>, normalVectors: ArrayLike<number>, binormalVectors: ArrayLike<number>, width: number, leftHeight: number, rightHeight: number) {
    const { vertices, normals, indices } = state;
    const vertexCount = vertices.elementCount;

    Vec3.fromArray(verticalLeftVector, normalVectors, offset);
    Vec3.scale(verticalLeftVector, verticalLeftVector, leftHeight);

    Vec3.fromArray(verticalRightVector, normalVectors, offset);
    Vec3.scale(verticalRightVector, verticalRightVector, rightHeight);

    Vec3.fromArray(horizontalVector, binormalVectors, offset);
    Vec3.scale(horizontalVector, horizontalVector, width);

    Vec3.fromArray(positionVector, controlPoints, offset);

    Vec3.add(p1, Vec3.add(p1, positionVector, horizontalVector), verticalRightVector);
    Vec3.sub(p2, Vec3.add(p2, positionVector, horizontalVector), verticalLeftVector);
    Vec3.sub(p3, Vec3.sub(p3, positionVector, horizontalVector), verticalLeftVector);
    Vec3.add(p4, Vec3.sub(p4, positionVector, horizontalVector), verticalRightVector);

    if (leftHeight < rightHeight) {
        ChunkedArray.add3(vertices, p4[0], p4[1], p4[2]);
        ChunkedArray.add3(vertices, p3[0], p3[1], p3[2]);
        ChunkedArray.add3(vertices, p2[0], p2[1], p2[2]);
        ChunkedArray.add3(vertices, p1[0], p1[1], p1[2]);
        Vec3.copy(verticalVector, verticalRightVector);
    } else {
        ChunkedArray.add3(vertices, p1[0], p1[1], p1[2]);
        ChunkedArray.add3(vertices, p2[0], p2[1], p2[2]);
        ChunkedArray.add3(vertices, p3[0], p3[1], p3[2]);
        ChunkedArray.add3(vertices, p4[0], p4[1], p4[2]);
        Vec3.copy(verticalVector, verticalLeftVector);
    }

    Vec3.cross(normalVector, horizontalVector, verticalVector);

    for (let i = 0; i < 4; ++i) {
        ChunkedArray.add3(normals, normalVector[0], normalVector[1], normalVector[2]);
    }
    ChunkedArray.add3(indices, vertexCount + 2, vertexCount + 1, vertexCount);
    ChunkedArray.add3(indices, vertexCount, vertexCount + 3, vertexCount + 2);
}

/** set arrowHeight = 0 for no arrow */
export function addSheet(state: MeshBuilder.State, controlPoints: ArrayLike<number>, normalVectors: ArrayLike<number>, binormalVectors: ArrayLike<number>, linearSegments: number, widthValues: ArrayLike<number>, heightValues: ArrayLike<number>, arrowHeight: number, startCap: boolean, endCap: boolean) {
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

        Vec3.add(tA, Vec3.add(tA, positionVector, horizontalVector), verticalVector);
        Vec3.copy(tB, normalVector);
        ChunkedArray.add3(vertices, tA[0], tA[1], tA[2]);
        ChunkedArray.add3(normals, tB[0], tB[1], tB[2]);

        Vec3.add(tA, Vec3.sub(tA, positionVector, horizontalVector), verticalVector);
        ChunkedArray.add3(vertices, tA[0], tA[1], tA[2]);
        ChunkedArray.add3(normals, tB[0], tB[1], tB[2]);

        // Vec3.add(tA, Vec3.sub(tA, positionVector, horizontalVector), verticalVector) // reuse tA
        Vec3.add(tB, Vec3.negate(tB, torsionVector), normalOffset);
        ChunkedArray.add3(vertices, tA[0], tA[1], tA[2]);
        ChunkedArray.add3(normals, tB[0], tB[1], tB[2]);

        Vec3.sub(tA, Vec3.sub(tA, positionVector, horizontalVector), verticalVector);
        ChunkedArray.add3(vertices, tA[0], tA[1], tA[2]);
        ChunkedArray.add3(normals, tB[0], tB[1], tB[2]);

        // Vec3.sub(tA, Vec3.sub(tA, positionVector, horizontalVector), verticalVector) // reuse tA
        Vec3.negate(tB, normalVector);
        ChunkedArray.add3(vertices, tA[0], tA[1], tA[2]);
        ChunkedArray.add3(normals, tB[0], tB[1], tB[2]);

        Vec3.sub(tA, Vec3.add(tA, positionVector, horizontalVector), verticalVector);
        ChunkedArray.add3(vertices, tA[0], tA[1], tA[2]);
        ChunkedArray.add3(normals, tB[0], tB[1], tB[2]);

        // Vec3.sub(tA, Vec3.add(tA, positionVector, horizontalVector), verticalVector) // reuse tA
        Vec3.add(tB, torsionVector, normalOffset);
        ChunkedArray.add3(vertices, tA[0], tA[1], tA[2]);
        ChunkedArray.add3(normals, tB[0], tB[1], tB[2]);

        Vec3.add(tA, Vec3.add(tA, positionVector, horizontalVector), verticalVector);
        ChunkedArray.add3(vertices, tA[0], tA[1], tA[2]);
        ChunkedArray.add3(normals, tB[0], tB[1], tB[2]);
    }

    for (let i = 0; i < linearSegments; ++i) {
        // the triangles are arranged such that opposing triangles of the sheet align
        // which prevents triangle intersection within tight curves
        for (let j = 0; j < 2; j++) {
            ChunkedArray.add3(
                indices,
                vertexCount + i * 8 + 2 * j, // a
                vertexCount + (i + 1) * 8 + 2 * j + 1, // c
                vertexCount + i * 8 + 2 * j + 1 // b
            );
            ChunkedArray.add3(
                indices,
                vertexCount + i * 8 + 2 * j, // a
                vertexCount + (i + 1) * 8 + 2 * j, // d
                vertexCount + (i + 1) * 8 + 2 * j + 1 // c
            );
        }
        for (let j = 2; j < 4; j++) {
            ChunkedArray.add3(
                indices,
                vertexCount + i * 8 + 2 * j, // a
                vertexCount + (i + 1) * 8 + 2 * j, // d
                vertexCount + i * 8 + 2 * j + 1, // b
            );
            ChunkedArray.add3(
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
        addCap(0, state, controlPoints, normalVectors, binormalVectors, width, h, h);
    } else if (arrowHeight > 0) {
        const width = widthValues[0];
        const height = heightValues[0];
        addCap(0, state, controlPoints, normalVectors, binormalVectors, width, arrowHeight, -height);
        addCap(0, state, controlPoints, normalVectors, binormalVectors, width, -arrowHeight, height);
    }

    if (endCap && arrowHeight === 0) {
        const width = widthValues[linearSegments];
        const height = heightValues[linearSegments];
        // use negative height to flip the direction the cap's triangles are facing
        addCap(linearSegments * 3, state, controlPoints, normalVectors, binormalVectors, width, -height, -height);
    }

    const addedVertexCount = (linearSegments + 1) * 8 +
        (startCap ? 4 : (arrowHeight > 0 ? 8 : 0)) +
        (endCap && arrowHeight === 0 ? 4 : 0);
    for (let i = 0, il = addedVertexCount; i < il; ++i) ChunkedArray.add(groups, currentGroup);
}