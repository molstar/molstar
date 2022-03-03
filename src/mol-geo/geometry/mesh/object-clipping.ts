/**
 * Copyright (c) 2021 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { ChunkedArray } from '../../../mol-data/util';
import { MeshValues } from '../../../mol-gl/renderable/mesh';
import { Sphere3D } from '../../../mol-math/geometry';
import { sphereIntersect } from '../../../mol-math/intersect';
import { Vec3 } from '../../../mol-math/linear-algebra';
import { ValueCell } from '../../../mol-util';
import { Clip } from '../../../mol-util/clip';

interface ObjectClippingInput {
    drawCount: number
    vertexCount: number
    instanceCount: number
    groupCount: number
    transformBuffer: Float32Array
    instanceBuffer: Float32Array
    positionBuffer: Float32Array
    indexBuffer: Uint32Array
    normalBuffer: Float32Array
    groupBuffer: Float32Array
    // colorData: TextureImage<Uint8Array>
    // colorType: 'group' | 'groupInstance'
    boundingSphere: Sphere3D
    invariantBoundingSphere: Sphere3D

    objectCount: number
    objectType: number[]
    objectPosition: number[]
    objectScale: number[]
    objectInvert: boolean[]
}

export function calcMeshObjectClipping(input: ObjectClippingInput) {
    console.log('calcMeshObjectClipping', input);
    const { drawCount, vertexCount, indexBuffer, positionBuffer, normalBuffer, objectCount, groupBuffer, objectPosition, objectScale, objectType, objectInvert } = input;
    const triangleCount = drawCount / 3;

    if (objectCount === 0) {
        return {
            indexBuffer,
            positionBuffer,
            normalBuffer,
            groupBuffer,
            drawCount,
            vertexCount
        };
    }

    const tests: ((p: Vec3) => boolean)[] = [];
    for (let i = 0; i < objectCount; ++i) {
        if (objectType[i] === Clip.Type.sphere) {
            const c = Vec3.fromArray(Vec3(), objectPosition, i * 3);
            const r = objectScale[i * 3] / 2;
            if (objectInvert[i]) {
                tests.push((p: Vec3) => Vec3.distance(p, c) <= r);
            } else {
                tests.push((p: Vec3) => Vec3.distance(p, c) >= r);
            }
        }
    }

    const invert = objectInvert[0];
    const center = Vec3.fromArray(Vec3(), objectPosition, 0);
    const radius = objectScale[0] / 2;

    // new
    const index = ChunkedArray.create(Uint32Array, 3, 1024, triangleCount);

    // re-use
    const vertex = ChunkedArray.create(Float32Array, 3, 1024, positionBuffer);
    vertex.currentIndex = vertexCount * 3;
    vertex.elementCount = vertexCount;
    const normal = ChunkedArray.create(Float32Array, 3, 1024, normalBuffer);
    normal.currentIndex = vertexCount * 3;
    normal.elementCount = vertexCount;
    const group = ChunkedArray.create(Float32Array, 1, 1024, groupBuffer);
    group.currentIndex = vertexCount;
    group.elementCount = vertexCount;

    const vA = Vec3();
    const vB = Vec3();
    const vC = Vec3();

    const xA = Vec3();
    const xB = Vec3();

    const nA = Vec3();
    const nB = Vec3();

    function add(a: number, b: number, p: Vec3) {
        Vec3.fromArray(nA, normalBuffer, a * 3);
        Vec3.fromArray(nB, normalBuffer, b * 3);
        Vec3.scale(nA, Vec3.add(nA, nA, nB), 0.5);
        ChunkedArray.add3(vertex, p[0], p[1], p[2]);
        ChunkedArray.add3(normal, nA[0], nA[1], nA[2]);
        ChunkedArray.add(group, groupBuffer[a]);
    }

    function split(tip: Vec3, sideA: Vec3, sideB: Vec3, tipIdx: number, sideIdxA: number, sideIdxB: number, reverse: boolean) {
        if (reverse) {
            sphereIntersect(xA, sideA, tip, center, radius, invert);
            sphereIntersect(xB, sideB, tip, center, radius, invert);
        } else {
            sphereIntersect(xA, tip, sideA, center, radius, invert);
            sphereIntersect(xB, tip, sideB, center, radius, invert);
        }
        add(tipIdx, sideIdxA, xA);
        add(tipIdx, sideIdxB, xB);
        if (reverse) {
            ChunkedArray.add3(index, tipIdx, newVertexCount + 1, newVertexCount);
            newTriangleCount += 1;
        } else {
            ChunkedArray.add3(index, newVertexCount, newVertexCount + 1, sideIdxA);
            ChunkedArray.add3(index, newVertexCount + 1, sideIdxB, sideIdxA);
            newTriangleCount += 2;
        }
        newVertexCount += 2;
    }

    let newVertexCount = vertexCount;
    let newTriangleCount = 0;

    for (let i = 0; i < triangleCount; ++i) {
        const iA = indexBuffer[i * 3];
        const iB = indexBuffer[i * 3 + 1];
        const iC = indexBuffer[i * 3 + 2];

        Vec3.fromArray(vA, positionBuffer, iA * 3);
        Vec3.fromArray(vB, positionBuffer, iB * 3);
        Vec3.fromArray(vC, positionBuffer, iC * 3);

        let tA = true;
        let tB = true;
        let tC = true;
        for (const t of tests) {
            if (!t(vA)) { tA = false; }
            if (!t(vB)) { tB = false; }
            if (!t(vC)) { tC = false; }
        }

        if (!tA && !tB && !tC) continue;

        if (tA && tB && tC) {
            ChunkedArray.add3(index, iA, iB, iC);
            newTriangleCount += 1;
        } else if (tA && tB && !tC) {
            split(vC, vA, vB, iC, iA, iB, false);
        } else if (!tA && !tB && tC) {
            split(vC, vA, vB, iC, iA, iB, true);
        } else if (tA && tC && !tB) {
            split(vB, vA, vC, iB, iA, iC, false);
        } else if (!tA && !tC && tB) {
            split(vB, vA, vC, iB, iA, iC, true);
        } else if (tB && tC && !tA) {
            split(vA, vC, vB, iA, iC, iB, false);
        } else if (!tB && !tC && tA) {
            split(vA, vC, vB, iA, iC, iB, true);
        } else {
            console.log('what');
        }
    }

    console.log({ vertexCount, newVertexCount, triangleCount, newTriangleCount });

    return {
        indexBuffer: ChunkedArray.compact(index),
        positionBuffer: ChunkedArray.compact(vertex),
        normalBuffer: ChunkedArray.compact(normal),
        groupBuffer: ChunkedArray.compact(group),
        drawCount: newTriangleCount * 3,
        vertexCount: newVertexCount
    };
}

//

export function applyMeshObjectClipping(values: MeshValues) {
    console.log('applyMeshObjectClipping');
    const clippingData = calcMeshObjectClipping({
        drawCount: values.drawCount.ref.value,
        vertexCount: values.uVertexCount.ref.value,
        instanceCount: values.instanceCount.ref.value,
        groupCount: values.uGroupCount.ref.value,
        transformBuffer: values.transform.ref.value,
        instanceBuffer: values.aInstance.ref.value,
        positionBuffer: values.aPosition.ref.value,
        indexBuffer: values.elements.ref.value,
        normalBuffer: values.aNormal.ref.value,
        groupBuffer: values.aGroup.ref.value,
        // colorData: TextureImage<Uint8Array>
        // colorType: 'group' | 'groupInstance'
        boundingSphere: values.boundingSphere.ref.value,
        invariantBoundingSphere: values.invariantBoundingSphere.ref.value,

        objectCount: values.dClipObjectCount.ref.value,
        objectType: values.uClipObjectType.ref.value,
        objectPosition: values.uClipObjectPosition.ref.value,
        objectScale: values.uClipObjectScale.ref.value,
        objectInvert: values.uClipObjectInvert.ref.value,
    });
    console.log(values.dClipObjectCount.ref.value, values.drawCount.ref.value, clippingData);

    ValueCell.updateIfChanged(values.dClipObjectCount, 0);
    ValueCell.update(values.elements, clippingData.indexBuffer);
    ValueCell.update(values.aPosition, clippingData.positionBuffer);
    ValueCell.update(values.aNormal, clippingData.normalBuffer);
    ValueCell.update(values.aGroup, clippingData.groupBuffer);
    ValueCell.updateIfChanged(values.drawCount, clippingData.drawCount);
    ValueCell.updateIfChanged(values.uVertexCount, clippingData.vertexCount);
}
