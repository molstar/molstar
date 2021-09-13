/**
 * Copyright (c) 2018-2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Vec3 } from '../../mol-math/linear-algebra';
import { Primitive, PrimitiveBuilder } from './primitive';
import { polygon } from './polygon';
import { Cage, createCage } from './cage';

const a = Vec3(), b = Vec3(), c = Vec3(), d = Vec3();
const points = polygon(4, true);

/**
 * Create a box
 */
function createBox(perforated: boolean): Primitive {
    const triangleCount = 12;
    const vertexCount = perforated ? 12 * 3 : 6 * 4;
    const builder = PrimitiveBuilder(triangleCount, vertexCount);

    // create sides
    for (let i = 0; i < 4; ++i) {
        const ni = (i + 1) % 4;
        Vec3.set(a, points[i * 3], points[i * 3 + 1], -0.5);
        Vec3.set(b, points[ni * 3], points[ni * 3 + 1], -0.5);
        Vec3.set(c, points[ni * 3], points[ni * 3 + 1], 0.5);
        Vec3.set(d, points[i * 3], points[i * 3 + 1], 0.5);
        if (perforated) {
            builder.add(a, b, c);
        } else {
            builder.addQuad(a, b, c, d);
        }
    }

    // create bases
    Vec3.set(a, points[0], points[1], -0.5);
    Vec3.set(b, points[3], points[4], -0.5);
    Vec3.set(c, points[6], points[7], -0.5);
    Vec3.set(d, points[9], points[10], -0.5);
    if (perforated) {
        builder.add(c, b, a);
    } else {
        builder.addQuad(d, c, b, a);
    }
    Vec3.set(a, points[0], points[1], 0.5);
    Vec3.set(b, points[3], points[4], 0.5);
    Vec3.set(c, points[6], points[7], 0.5);
    Vec3.set(d, points[9], points[10], 0.5);
    if (perforated) {
        builder.add(a, b, c);
    } else {
        builder.addQuad(a, b, c, d);
    }

    return builder.getPrimitive();
}

let box: Primitive;
export function Box() {
    if (!box) box = createBox(false);
    return box;
}

let perforatedBox: Primitive;
export function PerforatedBox() {
    if (!perforatedBox) perforatedBox = createBox(true);
    return perforatedBox;
}

let boxCage: Cage;
export function BoxCage() {
    if (!boxCage) {
        boxCage = createCage([
            0.5, 0.5, -0.5, // bottom
            -0.5, 0.5, -0.5,
            -0.5, -0.5, -0.5,
            0.5, -0.5, -0.5,
            0.5, 0.5, 0.5, // top
            -0.5, 0.5, 0.5,
            -0.5, -0.5, 0.5,
            0.5, -0.5, 0.5
        ], [
            0, 4, 1, 5, 2, 6, 3, 7, // sides
            0, 1, 1, 2, 2, 3, 3, 0, // bottom base
            4, 5, 5, 6, 6, 7, 7, 4 // top base
        ]);
    }
    return boxCage;
}