/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Vec3 } from 'mol-math/linear-algebra'
import { Primitive } from './primitive';

export const DefaultStarProps = {
    pointCount: 5,
    outerRadius: 1,
    innerRadius: 0.5,
    thickness: 0.3
}
export type StarProps = Partial<typeof DefaultStarProps>

const op = Vec3.zero()
const on = Vec3.zero()
const vn = Vec3.zero()
const p1 = Vec3.zero()
const p2 = Vec3.zero()
const p3 = Vec3.zero()

export function Star(props?: StarProps): Primitive {
    const { outerRadius, innerRadius, thickness, pointCount } = { ...DefaultStarProps, ...props }

    const triangleCount = pointCount * 2 * 2
    const vertexCount = triangleCount * 3

    const vertices = new Float32Array(vertexCount * 3)
    const normals = new Float32Array(vertexCount * 3)
    const indices = new Uint32Array(triangleCount * 3)

    const innerPoints = new Float32Array(pointCount * 2)
    const outerPoints = new Float32Array(pointCount * 2)

    for (let i = 0; i < pointCount; ++i) {
        const io = i * 2, ii = i * 2 + 1
        const co = io / pointCount * Math.PI, ci = ii / pointCount * Math.PI
        outerPoints[io] = Math.cos(co) * outerRadius
        outerPoints[ii] = Math.sin(co) * outerRadius
        innerPoints[io] = Math.cos(ci) * innerRadius
        innerPoints[ii] = Math.sin(ci) * innerRadius
    }

    Vec3.set(op, 0, 0, thickness / 2)
    Vec3.set(on, 0, 0, -thickness / 2)

    function add(a: Vec3, b: Vec3, c: Vec3, offset: number) {
        Vec3.toArray(a, vertices, offset)
        Vec3.toArray(b, vertices, offset + 3)
        Vec3.toArray(c, vertices, offset + 6)
        Vec3.triangleNormal(vn, a, b, c)
        for (let j = 0; j < 3; ++j) {
            Vec3.toArray(vn, normals, offset + 3 * j)
            indices[offset / 3 + j] = offset / 3 + j
        }
    }

    for (let i = 0; i < pointCount; ++i) {
        const ni = (i + 1) % pointCount
        Vec3.set(p1, outerPoints[i * 2], outerPoints[i * 2 + 1], 0)
        Vec3.set(p2, innerPoints[i * 2], innerPoints[i * 2 + 1], 0)
        Vec3.set(p3, outerPoints[ni * 2], outerPoints[ni * 2 + 1], 0)

        const offset = i * 3 * 3 * 4
        add(op, p1, p2, offset)
        add(on, p1, p2, offset + 9)
        add(op, p2, p3, offset + 18)
        add(on, p2, p3, offset + 27)
    }

    return { vertices, normals, indices }
}