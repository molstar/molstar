/**
 * Copyright (c) 2021 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

// adapted from three.js, MIT License Copyright 2010-2021 three.js authors

import { Vec3 } from '../../mol-math/linear-algebra';
import { Primitive } from './primitive';

export const DefaultTorusProps = {
    radius: 1,
    tube: 0.4,
    radialSegments: 8,
    tubularSegments: 6,
    arc: Math.PI * 2,
};
export type TorusProps = Partial<typeof DefaultTorusProps>

export function Torus(props?: TorusProps): Primitive {
    const { radius, tube, radialSegments, tubularSegments, arc } = { ...DefaultTorusProps, ...props };

    // buffers
    const indices: number[] = [];
    const vertices: number[] = [];
    const normals: number[] = [];

    // helper variables
    const center = Vec3();
    const vertex = Vec3();
    const normal = Vec3();

    // generate vertices and normals
    for (let j = 0; j <= radialSegments; ++j) {
        for (let i = 0; i <= tubularSegments; ++i) {
            const u = i / tubularSegments * arc;
            const v = j / radialSegments * Math.PI * 2;

            // vertex
            Vec3.set(
                vertex,
                (radius + tube * Math.cos(v)) * Math.cos(u),
                (radius + tube * Math.cos(v)) * Math.sin(u),
                tube * Math.sin(v)
            );
            vertices.push(...vertex);

            // normal
            Vec3.set(center, radius * Math.cos(u), radius * Math.sin(u), 0);
            Vec3.sub(normal, vertex, center);
            Vec3.normalize(normal, normal);
            normals.push(...normal);
        }
    }

    // generate indices
    for (let j = 1; j <= radialSegments; ++j) {
        for (let i = 1; i <= tubularSegments; ++i) {

            // indices
            const a = (tubularSegments + 1) * j + i - 1;
            const b = (tubularSegments + 1) * (j - 1) + i - 1;
            const c = (tubularSegments + 1) * (j - 1) + i;
            const d = (tubularSegments + 1) * j + i;

            // faces
            indices.push(a, b, d);
            indices.push(b, c, d);
        }
    }

    return {
        vertices: new Float32Array(vertices),
        normals: new Float32Array(normals),
        indices: new Uint32Array(indices)
    };
}