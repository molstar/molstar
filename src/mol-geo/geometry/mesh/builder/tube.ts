/**
 * Copyright (c) 2018-2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Vec3 } from '../../../../mol-math/linear-algebra';
import { ChunkedArray } from '../../../../mol-data/util';
import { MeshBuilder } from '../mesh-builder';

const normalVector = Vec3();
const surfacePoint = Vec3();
const controlPoint = Vec3();
const u = Vec3();
const v = Vec3();

function add2AndScale2(out: Vec3, a: Vec3, b: Vec3, sa: number, sb: number) {
    out[0] = (a[0] * sa) + (b[0] * sb);
    out[1] = (a[1] * sa) + (b[1] * sb);
    out[2] = (a[2] * sa) + (b[2] * sb);
}

function add3AndScale2(out: Vec3, a: Vec3, b: Vec3, c: Vec3, sa: number, sb: number) {
    out[0] = (a[0] * sa) + (b[0] * sb) + c[0];
    out[1] = (a[1] * sa) + (b[1] * sb) + c[1];
    out[2] = (a[2] * sa) + (b[2] * sb) + c[2];
}

// avoiding namespace lookup improved performance in Chrome (Aug 2020)
const v3fromArray = Vec3.fromArray;
const v3normalize = Vec3.normalize;
const v3negate = Vec3.negate;
const v3copy = Vec3.copy;
const v3cross = Vec3.cross;
const caAdd3 = ChunkedArray.add3;

const CosSinCache = new Map<number, { cos: number[], sin: number[] }>();
function getCosSin(radialSegments: number) {
    if (!CosSinCache.has(radialSegments)) {
        const cos: number[] = [];
        const sin: number[] = [];
        for (let j = 0; j < radialSegments; ++j) {
            const t = 2 * Math.PI * j / radialSegments;
            cos[j] = Math.cos(t);
            sin[j] = Math.sin(t);
        }
        CosSinCache.set(radialSegments, { cos, sin });
    }
    return CosSinCache.get(radialSegments)!;
}

export function addTube(state: MeshBuilder.State, controlPoints: ArrayLike<number>, normalVectors: ArrayLike<number>, binormalVectors: ArrayLike<number>, linearSegments: number, radialSegments: number, widthValues: ArrayLike<number>, heightValues: ArrayLike<number>, startCap: boolean, endCap: boolean) {
    const { currentGroup, vertices, normals, indices, groups } = state;

    let vertexCount = vertices.elementCount;

    const { cos, sin } = getCosSin(radialSegments);

    for (let i = 0; i <= linearSegments; ++i) {
        const i3 = i * 3;
        v3fromArray(u, normalVectors, i3);
        v3fromArray(v, binormalVectors, i3);
        v3fromArray(controlPoint, controlPoints, i3);

        const width = widthValues[i];
        const height = heightValues[i];

        for (let j = 0; j < radialSegments; ++j) {
            add3AndScale2(surfacePoint, u, v, controlPoint, height * cos[j], width * sin[j]);
            if (radialSegments === 2) {
                v3copy(normalVector, v);
                v3normalize(normalVector, normalVector);
                if (j !== 0 || i % 2 === 0) v3negate(normalVector, normalVector);
            } else {
                add2AndScale2(normalVector, u, v, width * cos[j], height * sin[j]);
            }
            v3normalize(normalVector, normalVector);

            caAdd3(vertices, surfacePoint[0], surfacePoint[1], surfacePoint[2]);
            caAdd3(normals, normalVector[0], normalVector[1], normalVector[2]);
        }
    }

    for (let i = 0; i < linearSegments; ++i) {
        for (let j = 0; j < radialSegments; ++j) {
            caAdd3(
                indices,
                vertexCount + i * radialSegments + (j + 1) % radialSegments,
                vertexCount + (i + 1) * radialSegments + (j + 1) % radialSegments,
                vertexCount + i * radialSegments + j
            );
            caAdd3(
                indices,
                vertexCount + (i + 1) * radialSegments + (j + 1) % radialSegments,
                vertexCount + (i + 1) * radialSegments + j,
                vertexCount + i * radialSegments + j
            );
        }
    }

    if (startCap) {
        const offset = 0;
        const centerVertex = vertices.elementCount;
        v3fromArray(u, normalVectors, offset);
        v3fromArray(v, binormalVectors, offset);
        v3fromArray(controlPoint, controlPoints, offset);
        v3cross(normalVector, v, u);

        caAdd3(vertices, controlPoint[0], controlPoint[1], controlPoint[2]);
        caAdd3(normals, normalVector[0], normalVector[1], normalVector[2]);

        const width = widthValues[0];
        const height = heightValues[0];

        vertexCount = vertices.elementCount;
        for (let i = 0; i < radialSegments; ++i) {
            add3AndScale2(surfacePoint, u, v, controlPoint, height * cos[i], width * sin[i]);

            caAdd3(vertices, surfacePoint[0], surfacePoint[1], surfacePoint[2]);
            caAdd3(normals, normalVector[0], normalVector[1], normalVector[2]);

            caAdd3(
                indices,
                vertexCount + (i + 1) % radialSegments,
                vertexCount + i,
                centerVertex
            );
        }
    }

    if (endCap) {
        const offset = linearSegments * 3;
        const centerVertex = vertices.elementCount;
        v3fromArray(u, normalVectors, offset);
        v3fromArray(v, binormalVectors, offset);
        v3fromArray(controlPoint, controlPoints, offset);
        v3cross(normalVector, u, v);

        caAdd3(vertices, controlPoint[0], controlPoint[1], controlPoint[2]);
        caAdd3(normals, normalVector[0], normalVector[1], normalVector[2]);

        const width = widthValues[linearSegments];
        const height = heightValues[linearSegments];

        vertexCount = vertices.elementCount;
        for (let i = 0; i < radialSegments; ++i) {
            add3AndScale2(surfacePoint, u, v, controlPoint, height * cos[i], width * sin[i]);

            caAdd3(vertices, surfacePoint[0], surfacePoint[1], surfacePoint[2]);
            caAdd3(normals, normalVector[0], normalVector[1], normalVector[2]);

            caAdd3(
                indices,
                vertexCount + i,
                vertexCount + (i + 1) % radialSegments,
                centerVertex
            );
        }
    }

    const addedVertexCount = (linearSegments + 1) * radialSegments + (startCap ? radialSegments + 1 : 0) + (endCap ? radialSegments + 1 : 0);
    ChunkedArray.addRepeat(groups, addedVertexCount, currentGroup);
}