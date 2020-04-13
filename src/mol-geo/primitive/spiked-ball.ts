/**
 * Copyright (c) 2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { createPrimitive, Primitive } from './primitive';
import { dodecahedronVertices, dodecahedronFaces } from './dodecahedron';
import { Vec3 } from '../../mol-math/linear-algebra';

function calcCenter(out: Vec3, ...vec3s: Vec3[]) {
    Vec3.set(out, 0, 0, 0);
    for (let i = 0, il = vec3s.length; i < il; ++i) {
        Vec3.add(out, out, vec3s[i]);
    }
    Vec3.scale(out, out, 1 / vec3s.length);
    return out;
}

const center = Vec3.zero();
const dir = Vec3.zero();
const tip = Vec3.zero();

const vecA = Vec3.zero();
const vecB = Vec3.zero();
const vecC = Vec3.zero();
const vecD = Vec3.zero();
const vecE = Vec3.zero();

/**
 * Create a spiked ball derived from a dodecahedron
 * @param radiusRatio ratio between inner radius (dodecahedron) and outher radius (spikes)
 */
export function SpikedBall(radiusRatio = 1): Primitive {
    const vertices = dodecahedronVertices.slice(0);
    const indices: number[] = [];

    let offset = vertices.length / 3;

    for (let i = 0, il = dodecahedronFaces.length; i < il; i += 5) {
        Vec3.fromArray(vecA, dodecahedronVertices, dodecahedronFaces[i] * 3);
        Vec3.fromArray(vecB, dodecahedronVertices, dodecahedronFaces[i + 1] * 3);
        Vec3.fromArray(vecC, dodecahedronVertices, dodecahedronFaces[i + 2] * 3);
        Vec3.fromArray(vecD, dodecahedronVertices, dodecahedronFaces[i + 3] * 3);
        Vec3.fromArray(vecE, dodecahedronVertices, dodecahedronFaces[i + 4] * 3);

        calcCenter(center, vecA, vecB, vecC, vecD, vecE);
        Vec3.triangleNormal(dir, vecA, vecB, vecC);
        Vec3.scaleAndAdd(tip, center, dir, radiusRatio);

        Vec3.toArray(tip, vertices, offset * 3);
        indices.push(offset, dodecahedronFaces[i], dodecahedronFaces[i + 1]);
        indices.push(offset, dodecahedronFaces[i + 1], dodecahedronFaces[i + 2]);
        indices.push(offset, dodecahedronFaces[i + 2], dodecahedronFaces[i + 3]);
        indices.push(offset, dodecahedronFaces[i + 3], dodecahedronFaces[i + 4]);
        indices.push(offset, dodecahedronFaces[i + 4], dodecahedronFaces[i]);

        offset += 1;
    }

    return createPrimitive(vertices, indices);
}