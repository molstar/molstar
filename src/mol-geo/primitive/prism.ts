/**
 * Copyright (c) 2018-2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Vec3 } from '../../mol-math/linear-algebra';
import { Primitive, PrimitiveBuilder } from './primitive';
import { polygon } from './polygon';
import { Cage } from './cage';

const on = Vec3(), op = Vec3();
const a = Vec3(), b = Vec3(), c = Vec3(), d = Vec3();

export const DefaultPrismProps = {
    height: 1,
    topCap: true,
    bottomCap: true,
};
export type PrismProps = Partial<typeof DefaultPrismProps>

/**
 * Create a prism with a base of 3 or more points
 */
export function Prism(points: ArrayLike<number>, props?: PrismProps): Primitive {
    const sideCount = points.length / 3;
    if (sideCount < 3) throw new Error('need at least 3 points to build a prism');

    const { height, topCap, bottomCap } = { ...DefaultPrismProps, ...props };

    let triangleCount = sideCount * 2;
    let vertexCount = sideCount * 4;

    const capCount = (topCap ? 1 : 0) + (bottomCap ? 1 : 0);
    if (sideCount === 3) {
        triangleCount += capCount;
        vertexCount += capCount * 3;
    } else if (sideCount === 4) {
        triangleCount += capCount * 2;
        vertexCount += capCount * 4;
    } else {
        triangleCount += capCount * sideCount;
        vertexCount += capCount * sideCount * 3;
    }

    const builder = PrimitiveBuilder(triangleCount, vertexCount);
    const halfHeight = height * 0.5;

    Vec3.set(on, 0, 0, -halfHeight);
    Vec3.set(op, 0, 0, halfHeight);

    // create sides
    for (let i = 0; i < sideCount; ++i) {
        const ni = (i + 1) % sideCount;
        Vec3.set(a, points[i * 3], points[i * 3 + 1], -halfHeight);
        Vec3.set(b, points[ni * 3], points[ni * 3 + 1], -halfHeight);
        Vec3.set(c, points[ni * 3], points[ni * 3 + 1], halfHeight);
        Vec3.set(d, points[i * 3], points[i * 3 + 1], halfHeight);
        builder.addQuad(a, b, c, d);
    }

    // create bases
    if (sideCount === 3) {
        if (topCap) {
            Vec3.set(a, points[0], points[1], -halfHeight);
            Vec3.set(b, points[3], points[4], -halfHeight);
            Vec3.set(c, points[6], points[7], -halfHeight);
            builder.add(c, b, a);
        }
        if (bottomCap) {
            Vec3.set(a, points[0], points[1], halfHeight);
            Vec3.set(b, points[3], points[4], halfHeight);
            Vec3.set(c, points[6], points[7], halfHeight);
            builder.add(a, b, c);
        }
    } else if (sideCount === 4) {
        if (topCap) {
            Vec3.set(a, points[0], points[1], -halfHeight);
            Vec3.set(b, points[3], points[4], -halfHeight);
            Vec3.set(c, points[6], points[7], -halfHeight);
            Vec3.set(d, points[9], points[10], -halfHeight);
            builder.addQuad(d, c, b, a);
        }
        if (bottomCap) {
            Vec3.set(a, points[0], points[1], halfHeight);
            Vec3.set(b, points[3], points[4], halfHeight);
            Vec3.set(c, points[6], points[7], halfHeight);
            Vec3.set(d, points[9], points[10], halfHeight);
            builder.addQuad(a, b, c, d);
        }
    } else {
        for (let i = 0; i < sideCount; ++i) {
            const ni = (i + 1) % sideCount;
            if (topCap) {
                Vec3.set(a, points[i * 3], points[i * 3 + 1], -halfHeight);
                Vec3.set(b, points[ni * 3], points[ni * 3 + 1], -halfHeight);
                builder.add(on, b, a);
            }
            if (bottomCap) {
                Vec3.set(a, points[i * 3], points[i * 3 + 1], halfHeight);
                Vec3.set(b, points[ni * 3], points[ni * 3 + 1], halfHeight);
                builder.add(a, b, op);
            }
        }
    }

    return builder.getPrimitive();
}

let diamond: Primitive;
export function DiamondPrism() {
    if (!diamond) diamond = Prism(polygon(4, false));
    return diamond;
}

let pentagonalPrism: Primitive;
export function PentagonalPrism() {
    if (!pentagonalPrism) pentagonalPrism = Prism(polygon(5, false));
    return pentagonalPrism;
}

let hexagonalPrism: Primitive;
export function HexagonalPrism() {
    if (!hexagonalPrism) hexagonalPrism = Prism(polygon(6, false));
    return hexagonalPrism;
}

let shiftedHexagonalPrism: Primitive;
export function ShiftedHexagonalPrism() {
    if (!shiftedHexagonalPrism) shiftedHexagonalPrism = Prism(polygon(6, true));
    return shiftedHexagonalPrism;
}

let heptagonalPrism: Primitive;
export function HeptagonalPrism() {
    if (!heptagonalPrism) heptagonalPrism = Prism(polygon(7, false));
    return heptagonalPrism;
}

//

/**
 * Create a prism cage
 */
export function PrismCage(points: ArrayLike<number>, height = 1): Cage {
    const sideCount = points.length / 3;

    const vertices: number[] = [];
    const edges: number[] = [];

    const halfHeight = height * 0.5;

    let offset = 0;

    // vertices and side edges
    for (let i = 0; i < sideCount; ++i) {
        vertices.push(
            points[i * 3], points[i * 3 + 1], -halfHeight,
            points[i * 3], points[i * 3 + 1], halfHeight
        );
        edges.push(offset, offset + 1);
        offset += 2;
    }

    // bases edges
    for (let i = 0; i < sideCount; ++i) {
        const ni = (i + 1) % sideCount;
        edges.push(
            i * 2, ni * 2,
            i * 2 + 1, ni * 2 + 1
        );
    }

    return { vertices, edges };
}

let diamondCage: Cage;
export function DiamondPrismCage() {
    if (!diamondCage) diamondCage = PrismCage(polygon(4, false));
    return diamondCage;
}

let pentagonalPrismCage: Cage;
export function PentagonalPrismCage() {
    if (!pentagonalPrismCage) pentagonalPrismCage = PrismCage(polygon(5, false));
    return pentagonalPrismCage;
}

let hexagonalPrismCage: Cage;
export function HexagonalPrismCage() {
    if (!hexagonalPrismCage) hexagonalPrismCage = PrismCage(polygon(6, false));
    return hexagonalPrismCage;
}