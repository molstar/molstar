/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Vec3 } from 'mol-math/linear-algebra'
import { Primitive, PrimitiveBuilder } from './primitive';
import { polygon } from './polygon'

const on = Vec3.create(0, 0, -0.5), op = Vec3.create(0, 0, 0.5)
const a = Vec3.zero(), b = Vec3.zero(), c = Vec3.zero(), d = Vec3.zero()

/**
 * Create a prism with a poligonal base
 */
export function Prism(points: ArrayLike<number>): Primitive {
    const sideCount = points.length / 2
    const baseCount = sideCount === 3 ? 1 : sideCount === 4 ? 2 : sideCount
    const count = 2 * baseCount + 2 * sideCount
    const builder = PrimitiveBuilder(count)

    // create sides
    for (let i = 0; i < sideCount; ++i) {
        const ni = (i + 1) % sideCount
        Vec3.set(a, points[i * 2], points[i * 2 + 1], -0.5)
        Vec3.set(b, points[ni * 2], points[ni * 2 + 1], -0.5)
        Vec3.set(c, points[ni * 2], points[ni * 2 + 1], 0.5)
        Vec3.set(d, points[i * 2], points[i * 2 + 1], 0.5)
        builder.add(a, b, c)
        builder.add(c, d, a)
    }

    // create bases
    if (sideCount === 3) {
        Vec3.set(a, points[0], points[1], -0.5)
        Vec3.set(b, points[2], points[3], -0.5)
        Vec3.set(c, points[4], points[5], -0.5)
        builder.add(a, b, c)
        Vec3.set(a, points[0], points[1], 0.5)
        Vec3.set(b, points[2], points[3], 0.5)
        Vec3.set(c, points[4], points[5], 0.5)
        builder.add(c, b, a)
    } else if (sideCount === 4) {
        Vec3.set(a, points[0], points[1], -0.5)
        Vec3.set(b, points[2], points[3], -0.5)
        Vec3.set(c, points[4], points[5], -0.5)
        Vec3.set(d, points[6], points[7], -0.5)
        builder.add(a, b, c)
        builder.add(c, d, a)
        Vec3.set(a, points[0], points[1], 0.5)
        Vec3.set(b, points[2], points[3], 0.5)
        Vec3.set(c, points[4], points[5], 0.5)
        Vec3.set(d, points[6], points[7], 0.5)
        builder.add(a, b, c)
        builder.add(c, d, a)
    } else {
        for (let i = 0; i < sideCount; ++i) {
            const ni = (i + 1) % sideCount
            Vec3.set(a, points[i * 2], points[i * 2 + 1], -0.5)
            Vec3.set(b, points[ni * 2], points[ni * 2 + 1], -0.5)
            builder.add(a, b, on)
            Vec3.set(a, points[i * 2], points[i * 2 + 1], 0.5)
            Vec3.set(b, points[ni * 2], points[ni * 2 + 1], 0.5)
            builder.add(op, b, a)
        }
    }

    return builder.getPrimitive()
}

let wedge: Primitive
export function Wedge() {
    if (!wedge) wedge = Prism(polygon(3, false))
    return wedge
}

let box: Primitive
export function Box() {
    if (!box) box = Prism(polygon(4, true))
    return box
}

let diamond: Primitive
export function DiamondPrism() {
    if (!diamond) diamond = Prism(polygon(4, false))
    return diamond
}

let pentagonalPrism: Primitive
export function PentagonalPrism() {
    if (!pentagonalPrism) pentagonalPrism = Prism(polygon(5, false))
    return pentagonalPrism
}

let hexagonalPrism: Primitive
export function HexagonalPrism() {
    if (!hexagonalPrism) hexagonalPrism = Prism(polygon(6, true))
    return hexagonalPrism
}