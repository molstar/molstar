/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Vec3 } from '../../mol-math/linear-algebra';
import { Primitive, PrimitiveBuilder } from './primitive';

export const DefaultStarProps = {
    pointCount: 5,
    outerRadius: 1,
    innerRadius: 0.5,
    thickness: 0.3
};
export type StarProps = Partial<typeof DefaultStarProps>

const op = Vec3.zero(), on = Vec3.zero();
const a = Vec3.zero(), b = Vec3.zero(), c = Vec3.zero();

export function Star(props?: StarProps): Primitive {
    const { outerRadius, innerRadius, thickness, pointCount } = { ...DefaultStarProps, ...props };

    const triangleCount = pointCount * 2 * 2;
    const builder = PrimitiveBuilder(triangleCount);

    const innerPoints = new Float32Array(pointCount * 2);
    const outerPoints = new Float32Array(pointCount * 2);

    for (let i = 0; i < pointCount; ++i) {
        const io = i * 2, ii = i * 2 + 1;
        const co = (io + 1) / pointCount * Math.PI, ci = (ii + 1) / pointCount * Math.PI;
        outerPoints[io] = Math.cos(co) * outerRadius;
        outerPoints[ii] = Math.sin(co) * outerRadius;
        innerPoints[io] = Math.cos(ci) * innerRadius;
        innerPoints[ii] = Math.sin(ci) * innerRadius;
    }

    Vec3.set(op, 0, 0, thickness / 2);
    Vec3.set(on, 0, 0, -thickness / 2);

    for (let i = 0; i < pointCount; ++i) {
        const ni = (i + 1) % pointCount;
        Vec3.set(a, outerPoints[i * 2], outerPoints[i * 2 + 1], 0);
        Vec3.set(b, innerPoints[i * 2], innerPoints[i * 2 + 1], 0);
        Vec3.set(c, outerPoints[ni * 2], outerPoints[ni * 2 + 1], 0);

        builder.add(op, a, b);
        builder.add(b, a, on);
        builder.add(op, b, c);
        builder.add(c, b, on);
    }

    return builder.getPrimitive();
}