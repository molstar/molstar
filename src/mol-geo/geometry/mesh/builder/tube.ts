/**
 * Copyright (c) 2018-2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Vec3 } from '../../../../mol-math/linear-algebra';
import { ChunkedArray } from '../../../../mol-data/util';
import { MeshBuilder } from '../mesh-builder';

const normalVector = Vec3.zero();
const surfacePoint = Vec3.zero();
const controlPoint = Vec3.zero();
const u = Vec3.zero();
const v = Vec3.zero();

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

export function addTube(state: MeshBuilder.State, controlPoints: ArrayLike<number>, normalVectors: ArrayLike<number>, binormalVectors: ArrayLike<number>, linearSegments: number, radialSegments: number, widthValues: ArrayLike<number>, heightValues: ArrayLike<number>, waveFactor: number, startCap: boolean, endCap: boolean) {
    const { currentGroup, vertices, normals, indices, groups } = state;

    let vertexCount = vertices.elementCount;
    const di = 1 / linearSegments;

    for (let i = 0; i <= linearSegments; ++i) {
        const i3 = i * 3;
        Vec3.fromArray(u, normalVectors, i3);
        Vec3.fromArray(v, binormalVectors, i3);
        Vec3.fromArray(controlPoint, controlPoints, i3);

        const width = widthValues[i];
        const height = heightValues[i];

        const tt = di * i - 0.5;
        const ff = 1 + (waveFactor - 1) * (Math.cos(2 * Math.PI * tt) + 1);
        const w = ff * width, h = ff * height;

        for (let j = 0; j < radialSegments; ++j) {
            const t = 2 * Math.PI * j / radialSegments;

            add3AndScale2(surfacePoint, u, v, controlPoint, h * Math.cos(t), w * Math.sin(t));
            if (radialSegments === 2) {
                // add2AndScale2(normalVector, u, v, w * Math.cos(t), h * Math.sin(t))
                Vec3.copy(normalVector, v);
                Vec3.normalize(normalVector, normalVector);
                if (t !== 0 || i % 2 === 0) Vec3.negate(normalVector, normalVector);
            } else {
                add2AndScale2(normalVector, u, v, w * Math.cos(t), h * Math.sin(t));
            }
            Vec3.normalize(normalVector, normalVector);

            ChunkedArray.add3(vertices, surfacePoint[0], surfacePoint[1], surfacePoint[2]);
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
        const offset = 0;
        const centerVertex = vertices.elementCount;
        Vec3.fromArray(u, normalVectors, offset);
        Vec3.fromArray(v, binormalVectors, offset);
        Vec3.fromArray(controlPoint, controlPoints, offset);
        Vec3.cross(normalVector, v, u);

        ChunkedArray.add3(vertices, controlPoint[0], controlPoint[1], controlPoint[2]);
        ChunkedArray.add3(normals, normalVector[0], normalVector[1], normalVector[2]);

        const width = widthValues[0];
        const height = heightValues[0];

        vertexCount = vertices.elementCount;
        for (let i = 0; i < radialSegments; ++i) {
            const t = 2 * Math.PI * i / radialSegments;

            add3AndScale2(surfacePoint, u, v, controlPoint, height * Math.cos(t), width * Math.sin(t));

            ChunkedArray.add3(vertices, surfacePoint[0], surfacePoint[1], surfacePoint[2]);
            ChunkedArray.add3(normals, normalVector[0], normalVector[1], normalVector[2]);

            ChunkedArray.add3(
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
        Vec3.fromArray(u, normalVectors, offset);
        Vec3.fromArray(v, binormalVectors, offset);
        Vec3.fromArray(controlPoint, controlPoints, offset);
        Vec3.cross(normalVector, u, v);

        ChunkedArray.add3(vertices, controlPoint[0], controlPoint[1], controlPoint[2]);
        ChunkedArray.add3(normals, normalVector[0], normalVector[1], normalVector[2]);

        const width = widthValues[linearSegments];
        const height = heightValues[linearSegments];

        vertexCount = vertices.elementCount;
        for (let i = 0; i < radialSegments; ++i) {
            const t = 2 * Math.PI * i / radialSegments;

            add3AndScale2(surfacePoint, u, v, controlPoint, height * Math.cos(t), width * Math.sin(t));

            ChunkedArray.add3(vertices, surfacePoint[0], surfacePoint[1], surfacePoint[2]);
            ChunkedArray.add3(normals, normalVector[0], normalVector[1], normalVector[2]);

            ChunkedArray.add3(
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