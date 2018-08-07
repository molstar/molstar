/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Vec3 } from 'mol-math/linear-algebra'
import { Primitive } from './primitive';

export const DefaultWedgeProps = {
    width: 1,
    height: 1,
    depth: 1
}
export type WedgeProps = Partial<typeof DefaultWedgeProps>

const _a = Vec3.create(0, 0.5, 0.5)
const _b = Vec3.create(0.5, -0.5, 0.5)
const _c = Vec3.create(-0.5, -0.5, 0.5)
const _d = Vec3.create(0, 0.5, -0.5)
const _e = Vec3.create(0.5, -0.5, -0.5)
const _f = Vec3.create(-0.5, -0.5, -0.5)

const a = Vec3.zero(), b = Vec3.zero(), c = Vec3.zero()
const d = Vec3.zero(), e = Vec3.zero(), f = Vec3.zero()

const nabc = Vec3.create(0, 0, 1)
const ndef = Vec3.create(0, 0, -1)
const nabde = Vec3.zero()
const nbcef = Vec3.create(0, -1, 0)
const nacdf = Vec3.zero()

const s = Vec3.zero()

export function Wedge(props?: WedgeProps): Primitive {
    const { width, height, depth } = { ...DefaultWedgeProps, ...props }

    const vertices = new Float32Array(54)
    const normals = new Float32Array(54)
    const indices = new Uint32Array(24)

    Vec3.set(s, width, height, depth)
    Vec3.mul(a, _a, s); Vec3.mul(b, _b, s); Vec3.mul(c, _c, s)
    Vec3.mul(d, _d, s); Vec3.mul(e, _e, s); Vec3.mul(f, _f, s)

    Vec3.sub(nabde, b, a)
    Vec3.normalize(nabde, Vec3.set(nabde, -nabde[1], nabde[0], 0))
    Vec3.sub(nacdf, c, a)
    Vec3.normalize(nacdf, Vec3.set(nacdf, nacdf[1], -nacdf[0], 0))

    let vc = 0
    let ic = 0

    // abc
    Vec3.toArray(a, vertices, vc + 0)
    Vec3.toArray(c, vertices, vc + 3)
    Vec3.toArray(b, vertices, vc + 6)
    for (let i = 0; i < 3; ++i) Vec3.toArray(nabc, normals, vc + i * 3)
    indices[ic + 0] = vc / 3 + 0
    indices[ic + 1] = vc / 3 + 1
    indices[ic + 2] = vc / 3 + 2
    vc += 9
    ic += 3

    // def
    Vec3.toArray(d, vertices, vc + 0)
    Vec3.toArray(e, vertices, vc + 3)
    Vec3.toArray(f, vertices, vc + 6)
    for (let i = 0; i < 3; ++i) Vec3.toArray(ndef, normals, vc + i * 3)
    indices[ic + 0] = vc / 3 + 0
    indices[ic + 1] = vc / 3 + 1
    indices[ic + 2] = vc / 3 + 2
    vc += 9
    ic += 3

    // abde
    Vec3.toArray(a, vertices, vc + 0)
    Vec3.toArray(d, vertices, vc + 3)
    Vec3.toArray(e, vertices, vc + 6)
    Vec3.toArray(b, vertices, vc + 9)
    for (let i = 0; i < 4; ++i) Vec3.toArray(nabde, normals, vc + i * 3)
    indices[ic + 0] = vc / 3 + 2
    indices[ic + 1] = vc / 3 + 1
    indices[ic + 2] = vc / 3 + 0
    indices[ic + 3] = vc / 3 + 0
    indices[ic + 4] = vc / 3 + 3
    indices[ic + 5] = vc / 3 + 2
    vc += 12
    ic += 6

    // acdf
    Vec3.toArray(d, vertices, vc + 0)
    Vec3.toArray(a, vertices, vc + 3)
    Vec3.toArray(c, vertices, vc + 6)
    Vec3.toArray(f, vertices, vc + 9)
    for (let i = 0; i < 4; ++i) Vec3.toArray(nacdf, normals, vc + i * 3)
    indices[ic + 0] = vc / 3 + 2
    indices[ic + 1] = vc / 3 + 1
    indices[ic + 2] = vc / 3 + 0
    indices[ic + 3] = vc / 3 + 0
    indices[ic + 4] = vc / 3 + 3
    indices[ic + 5] = vc / 3 + 2
    vc += 12
    ic += 6

    // bcef
    Vec3.toArray(e, vertices, vc + 0)
    Vec3.toArray(f, vertices, vc + 3)
    Vec3.toArray(c, vertices, vc + 6)
    Vec3.toArray(b, vertices, vc + 9)
    for (let i = 0; i < 4; ++i) Vec3.toArray(nbcef, normals, vc + i * 3)
    indices[ic + 0] = vc / 3 + 2
    indices[ic + 1] = vc / 3 + 1
    indices[ic + 2] = vc / 3 + 0
    indices[ic + 3] = vc / 3 + 0
    indices[ic + 4] = vc / 3 + 3
    indices[ic + 5] = vc / 3 + 2
    vc += 12
    ic += 6

    return { vertices, normals, indices }
}