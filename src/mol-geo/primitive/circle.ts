/**
 * Copyright (c) 2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Primitive } from './primitive';
import { Cage } from './cage';

export const DefaultCircleProps = {
    radius: 1,
    segments: 4,
    thetaStart: 0,
    thetaLength: Math.PI * 2
}
export type CirclerProps = Partial<typeof DefaultCircleProps>

export function Circle(props?: CirclerProps): Primitive {
    const { radius, segments, thetaStart, thetaLength } = { ...DefaultCircleProps, ...props }

    const isFull = thetaLength === Math.PI * 2
    const count = isFull ? segments + 1 : segments + 2

    const vertices = new Float32Array(count * 3)
    const normals = new Float32Array(count * 3)
    const indices = new Uint32Array(segments * 3)

    // center
    vertices[0] = 0; vertices[1] = 0; vertices[2] = 0;
    normals[0] = 0; normals[1] = 0; normals[2] = 1;

    // vertices & normals
    for (let s = 0, i = 3; s < segments; ++s, i += 3) {
        const segment = thetaStart + s / segments * thetaLength;

        vertices[i] = radius * Math.cos(segment)
        vertices[i + 1] = radius * Math.sin(segment)
        vertices[i + 2] = 0

        normals[i] = 0; normals[i + 1] = 0; normals[i + 2] = 1;
    }

    // indices
    for (let s = 1, i = 0; s <= segments; ++s, i += 3) {
        indices[i] = s; indices[i + 1] = s + 1; indices[i + 2] = 0;
    }

    if (isFull) {
        // indices[i] = s; indices[i + 1] = s + 1; indices[i + 2] = 0;
    } else {
        const segment = thetaStart + thetaLength;
        const i = (segments + 1) * 3

        vertices[i] = radius * Math.cos(segment)
        vertices[i + 1] = radius * Math.sin(segment)
        vertices[i + 2] = 0

        normals[i] = 0; normals[i + 1] = 0; normals[i + 2] = 1;

        const j = segments * 3
        indices[j] = segments + 1
        indices[j + 1] = segments + 2
        indices[j + 2] = 0
    }

    return { vertices, normals, indices }
}

export function CircleCage(props?: CirclerProps): Cage {
    const { radius, segments, thetaStart, thetaLength } = { ...DefaultCircleProps, ...props }

    const n = segments * 3
    const vertices = new Float32Array(n)
    const edges = new Uint32Array(n)

    for (let s = 0, i = 0; s <= segments; ++s, i += 3) {
        const segment = thetaStart + s / segments * thetaLength;

        // vertex
        vertices[i] = radius * Math.cos(segment)
        vertices[i + 1] = radius * Math.sin(segment)
        vertices[i + 2] = 0
    }

    // indices
    for (let s = 1, i = 0; s <= segments; ++s, i += 3) {
        edges[0] = s; edges[1] = s + 1;
    }

    return { vertices, edges }
}