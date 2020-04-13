/**
 * Copyright (c) 2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Vec3, Mat4 } from '../../../../mol-math/linear-algebra';
import { MeshBuilder } from '../mesh-builder';
import { Axes3D } from '../../../../mol-math/geometry';
import { createCage } from '../../../primitive/cage';

const tmpVec = Vec3();
const tmpMatrix = Mat4.identity();

const tmpVertices = new Float32Array(6 * 3);
const tmpEdges = new Uint8Array([0, 1, 2, 3, 4, 5]);

export function addAxes(state: MeshBuilder.State, axes: Axes3D, radiusScale: number, detail: number, radialSegments: number) {
    const { origin, dirA, dirB, dirC } = axes;

    Vec3.add(tmpVec, origin, dirA);
    Vec3.toArray(Vec3.add(tmpVec, origin, dirA), tmpVertices, 0);
    Vec3.toArray(Vec3.sub(tmpVec, origin, dirA), tmpVertices, 3);
    Vec3.toArray(Vec3.add(tmpVec, origin, dirB), tmpVertices, 6);
    Vec3.toArray(Vec3.sub(tmpVec, origin, dirB), tmpVertices, 9);
    Vec3.toArray(Vec3.add(tmpVec, origin, dirC), tmpVertices, 12);
    Vec3.toArray(Vec3.sub(tmpVec, origin, dirC), tmpVertices, 15);

    const cage = createCage(tmpVertices, tmpEdges);
    const volume = Axes3D.volume(axes);
    const radius = (Math.cbrt(volume) / 300) * radiusScale;

    MeshBuilder.addCage(state, tmpMatrix, cage, radius, detail, radialSegments);
}