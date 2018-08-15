/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Vec3 } from 'mol-math/linear-algebra'
import { Primitive, PrimitiveBuilder } from './primitive';
import { polygon } from './polygon'

const a = Vec3.zero(), b = Vec3.zero(), c = Vec3.zero(), d = Vec3.zero()
const points = polygon(4, true)

/**
 * Create a box
 */
function createBox(perforated: boolean): Primitive {
    const builder = PrimitiveBuilder(12)

    // create sides
    for (let i = 0; i < 4; ++i) {
        const ni = (i + 1) % 4
        Vec3.set(a, points[i * 2], points[i * 2 + 1], -0.5)
        Vec3.set(b, points[ni * 2], points[ni * 2 + 1], -0.5)
        Vec3.set(c, points[ni * 2], points[ni * 2 + 1], 0.5)
        Vec3.set(d, points[i * 2], points[i * 2 + 1], 0.5)
        builder.add(a, b, c)
        if (!perforated) builder.add(c, d, a)
    }

    // create bases
    Vec3.set(a, points[0], points[1], -0.5)
    Vec3.set(b, points[2], points[3], -0.5)
    Vec3.set(c, points[4], points[5], -0.5)
    Vec3.set(d, points[6], points[7], -0.5)
    builder.add(a, b, c)
    if (!perforated) builder.add(c, d, a)
    Vec3.set(a, points[0], points[1], 0.5)
    Vec3.set(b, points[2], points[3], 0.5)
    Vec3.set(c, points[4], points[5], 0.5)
    Vec3.set(d, points[6], points[7], 0.5)
    builder.add(a, b, c)
    if (!perforated) builder.add(c, d, a)

    return builder.getPrimitive()
}

let box: Primitive
export function Box() {
    if (!box) box = createBox(false)
    return box
}

let perforatedBox: Primitive
export function PerforatedBox() {
    if (!perforatedBox) perforatedBox = createBox(true)
    return perforatedBox
}