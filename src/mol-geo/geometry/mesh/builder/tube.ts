/**
 * Copyright (c) 2018-2024 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Gianluca Tomasello <giagitom@gmail.com>
 */

import { Vec3 } from '../../../../mol-math/linear-algebra';
import { cantorPairing, ChunkedArray } from '../../../../mol-data/util';
import { MeshBuilder } from '../mesh-builder';

const normalVector = Vec3();
const capNormalSmoothingVector = Vec3();
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
const v3scaleAndAdd = Vec3.scaleAndAdd;
const v3cross = Vec3.cross;
const v3slerp = Vec3.slerp;
const v3dot = Vec3.dot;
const v3unitX = Vec3.unitX;
const caAdd3 = ChunkedArray.add3;

const CosSinCache = new Map<number, { cos: number[], sin: number[] }>();
function getCosSin(radialSegments: number, shift: boolean) {
    const offset = shift ? 1 : 0;
    const hash = cantorPairing(radialSegments, offset);
    if (!CosSinCache.has(hash)) {
        const cos: number[] = [];
        const sin: number[] = [];
        for (let j = 0; j < radialSegments; ++j) {
            const t = (j * 2 + offset) / radialSegments * Math.PI;
            cos[j] = Math.cos(t);
            sin[j] = Math.sin(t);
        }
        CosSinCache.set(hash, { cos, sin });
    }
    return CosSinCache.get(hash)!;
}

export function addTube(state: MeshBuilder.State, controlPoints: ArrayLike<number>, normalVectors: ArrayLike<number>, binormalVectors: ArrayLike<number>, linearSegments: number, radialSegments: number, widthValues: ArrayLike<number>, heightValues: ArrayLike<number>, startCap: boolean, endCap: boolean, crossSection: 'elliptical' | 'rounded', roundCap = false) {
    const { currentGroup, vertices, normals, indices, groups } = state;

    let vertexCount = vertices.elementCount;

    const { cos, sin } = getCosSin(radialSegments, crossSection === 'rounded');

    const q1 = Math.round(radialSegments / 4);
    const q3 = q1 * 3;

    const roundCapFlag = roundCap && linearSegments && (startCap || endCap);
    let halfLinearSegments;
    const doubleRoundCap = roundCapFlag && startCap && endCap;
    if (doubleRoundCap) halfLinearSegments = linearSegments / 2;
    for (let i = 0; i <= linearSegments; ++i) {
        const i3 = i * 3;
        v3fromArray(u, normalVectors, i3);
        v3fromArray(v, binormalVectors, i3);
        v3fromArray(controlPoint, controlPoints, i3);

        let width = widthValues[i];
        let height = heightValues[i];
        let capSmoothingFactor: number;
        if (roundCapFlag) {
            const sc = doubleRoundCap ? i <= halfLinearSegments! : startCap;
            if (doubleRoundCap) {
                capSmoothingFactor = Math.max(Number.EPSILON, Math.sqrt(1 - Math.pow((sc ? halfLinearSegments! - i : i - halfLinearSegments!) / halfLinearSegments!, 2)));
            } else {
                capSmoothingFactor = Math.max(Number.EPSILON, Math.sqrt(1 - Math.pow((sc ? linearSegments - i : i) / linearSegments, 2)));
            }
            width *= capSmoothingFactor;
            height *= capSmoothingFactor;
            v3cross(capNormalSmoothingVector, sc ? v : u, sc ? u : v);
            v3normalize(capNormalSmoothingVector, capNormalSmoothingVector);
        }
        const rounded = crossSection === 'rounded' && height > width;

        for (let j = 0; j < radialSegments; ++j) {
            if (rounded) {
                add3AndScale2(surfacePoint, u, v, controlPoint, width * cos[j], width * sin[j]);
                const h = v3dot(v, v3unitX) < 0
                    ? (j < q1 || j >= q3) ? height - width : -height + width
                    : (j >= q1 && j < q3) ? -height + width : height - width;
                v3scaleAndAdd(surfacePoint, surfacePoint, u, h);
                if (j === q1 || j === q1 - 1) {
                    add2AndScale2(normalVector, u, v, 0, 1);
                } else if (j === q3 || j === q3 - 1) {
                    add2AndScale2(normalVector, u, v, 0, -1);
                } else {
                    add2AndScale2(normalVector, u, v, cos[j], sin[j]);
                }
            } else {
                add3AndScale2(surfacePoint, u, v, controlPoint, height * cos[j], width * sin[j]);
                add2AndScale2(normalVector, u, v, width * cos[j], height * sin[j]);
            }
            v3normalize(normalVector, normalVector);

            caAdd3(vertices, surfacePoint[0], surfacePoint[1], surfacePoint[2]);
            if (roundCapFlag) {
                v3slerp(normalVector, capNormalSmoothingVector, normalVector, capSmoothingFactor!);
            }
            caAdd3(normals, normalVector[0], normalVector[1], normalVector[2]);
        }
    }

    const radialSegmentsHalf = Math.round(radialSegments / 2);

    for (let i = 0; i < linearSegments; ++i) {
        // the triangles are arranged such that opposing triangles of the sheet align
        // which prevents triangle intersection within tight curves
        for (let j = 0; j < radialSegmentsHalf; ++j) {
            caAdd3(
                indices,
                vertexCount + i * radialSegments + (j + 1) % radialSegments, // a
                vertexCount + (i + 1) * radialSegments + (j + 1) % radialSegments, // c
                vertexCount + i * radialSegments + j // b
            );
            caAdd3(
                indices,
                vertexCount + (i + 1) * radialSegments + (j + 1) % radialSegments, // c
                vertexCount + (i + 1) * radialSegments + j, // d
                vertexCount + i * radialSegments + j // b
            );
        }
        for (let j = radialSegmentsHalf; j < radialSegments; ++j) {
            caAdd3(
                indices,
                vertexCount + i * radialSegments + (j + 1) % radialSegments, // a
                vertexCount + (i + 1) * radialSegments + j, // d
                vertexCount + i * radialSegments + j // b
            );
            caAdd3(
                indices,
                vertexCount + (i + 1) * radialSegments + (j + 1) % radialSegments, // c
                vertexCount + (i + 1) * radialSegments + j, // d
                vertexCount + i * radialSegments + (j + 1) % radialSegments, // a
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

        const width = roundCapFlag ? 0 : widthValues[0];
        let height = roundCapFlag ? 0 : heightValues[0];
        const rounded = crossSection === 'rounded' && height > width;
        if (rounded) height -= width;

        vertexCount = vertices.elementCount;
        for (let i = 0; i < radialSegments; ++i) {
            if (rounded) {
                add3AndScale2(surfacePoint, u, v, controlPoint, width * cos[i], width * sin[i]);
                v3scaleAndAdd(surfacePoint, surfacePoint, u, (i < q1 || i >= q3) ? height : -height);
            } else {
                add3AndScale2(surfacePoint, u, v, controlPoint, height * cos[i], width * sin[i]);
            }

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

        const width = roundCapFlag ? 0 : widthValues[linearSegments];
        let height = roundCapFlag ? 0 : heightValues[linearSegments];
        const rounded = crossSection === 'rounded' && height > width;
        if (rounded) height -= width;

        vertexCount = vertices.elementCount;
        for (let i = 0; i < radialSegments; ++i) {
            if (rounded) {
                add3AndScale2(surfacePoint, u, v, controlPoint, width * cos[i], width * sin[i]);
                v3scaleAndAdd(surfacePoint, surfacePoint, u, (i < q1 || i >= q3) ? height : -height);
            } else {
                add3AndScale2(surfacePoint, u, v, controlPoint, height * cos[i], width * sin[i]);
            }

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