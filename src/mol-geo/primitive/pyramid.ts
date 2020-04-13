/**
 * Copyright (c) 2018-2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Vec3 } from '../../mol-math/linear-algebra';
import { Primitive, PrimitiveBuilder, createPrimitive } from './primitive';
import { polygon } from './polygon';
import { Cage } from './cage';

const on = Vec3.create(0, 0, -0.5), op = Vec3.create(0, 0, 0.5);
const a = Vec3.zero(), b = Vec3.zero(), c = Vec3.zero(), d = Vec3.zero();

/**
 * Create a pyramid with a polygonal base
 */
export function Pyramid(points: ArrayLike<number>): Primitive {
    const sideCount = points.length / 3;
    const baseCount = sideCount === 3 ? 1 : sideCount === 4 ? 2 : sideCount;
    const count = 2 * baseCount + 2 * sideCount;
    const builder = PrimitiveBuilder(count);

    // create sides
    for (let i = 0; i < sideCount; ++i) {
        const ni = (i + 1) % sideCount;
        Vec3.set(a, points[i * 3], points[i * 3 + 1], -0.5);
        Vec3.set(b, points[ni * 3], points[ni * 3 + 1], -0.5);
        builder.add(a, b, op);
    }

    // create base
    if (sideCount === 3) {
        Vec3.set(a, points[0], points[1], -0.5);
        Vec3.set(b, points[3], points[4], -0.5);
        Vec3.set(c, points[6], points[7], -0.5);
        builder.add(c, b, a);
    } else if (sideCount === 4) {
        Vec3.set(a, points[0], points[1], -0.5);
        Vec3.set(b, points[3], points[4], -0.5);
        Vec3.set(c, points[6], points[7], -0.5);
        Vec3.set(d, points[9], points[10], -0.5);
        builder.add(c, b, a);
        builder.add(a, d, c);
    } else {
        for (let i = 0; i < sideCount; ++i) {
            const ni = (i + 1) % sideCount;
            Vec3.set(a, points[i * 3], points[i * 3 + 1], -0.5);
            Vec3.set(b, points[ni * 3], points[ni * 3 + 1], -0.5);
            builder.add(on, b, a);
        }
    }

    return builder.getPrimitive();
}

let triangularPyramid: Primitive;
export function TriangularPyramid() {
    if (!triangularPyramid) triangularPyramid = Pyramid(polygon(3, true));
    return triangularPyramid;
}

let octagonalPyramid: Primitive;
export function OctagonalPyramid() {
    if (!octagonalPyramid) octagonalPyramid = Pyramid(polygon(8, true));
    return octagonalPyramid;
}

let perforatedOctagonalPyramid: Primitive;
export function PerforatedOctagonalPyramid() {
    if (!perforatedOctagonalPyramid) {
        const points = polygon(8, true);
        const vertices = new Float32Array(8 * 3 + 6);
        for (let i = 0; i < 8; ++i) {
            vertices[i * 3] = points[i * 3];
            vertices[i * 3 + 1] = points[i * 3 + 1];
            vertices[i * 3 + 2] = -0.5;
        }
        vertices[8 * 3] = 0;
        vertices[8 * 3 + 1] = 0;
        vertices[8 * 3 + 2] = -0.5;
        vertices[8 * 3 + 3] = 0;
        vertices[8 * 3 + 4] = 0;
        vertices[8 * 3 + 5] = 0.5;
        const indices: ReadonlyArray<number> = [
            0, 1, 8,  1, 2, 8,  4, 5, 8,  5, 6, 8,
            2, 3, 9,  3, 4, 9,  6, 7, 9,  7, 0, 9
        ];
        perforatedOctagonalPyramid = createPrimitive(vertices, indices);
    }
    return perforatedOctagonalPyramid;
}

//

/**
 * Create a prism cage
 */
export function PyramidCage(points: ArrayLike<number>): Cage {
    const sideCount = points.length / 3;

    // const count = 4 * sideCount
    const vertices: number[] = [];
    const edges: number[] = [];

    let offset = 1;
    vertices.push(op[0], op[1], op[2]);

    // vertices and side edges
    for (let i = 0; i < sideCount; ++i) {
        vertices.push(points[i * 3], points[i * 3 + 1], -0.5);
        edges.push(0, offset);
        offset += 1;
    }

    // bases edges
    for (let i = 0; i < sideCount; ++i) {
        const ni = (i + 1) % sideCount;
        edges.push(i + 1, ni + 1);
    }

    return { vertices, edges };
}

let octagonalPyramidCage: Cage;
export function OctagonalPyramidCage() {
    if (!octagonalPyramidCage) octagonalPyramidCage = PyramidCage(polygon(8, true));
    return octagonalPyramidCage;
}