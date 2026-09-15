/**
 * Copyright (c) 2018-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author Adam Midlik <midlik@gmail.com>
 */

import { Icosahedron } from './icosahedron';
import { Polyhedron } from './polyhedron';
import { Primitive } from './primitive';

/** Calculate vertex count for subdived icosahedron */
export function sphereVertexCount(detail: number) {
    return 10 * Math.pow(Math.pow(2, detail), 2) + 2;
}

/** Create sphere by subdividing an icosahedron */
export function Sphere(detail: number, options?: { subset: 'ring' | 'caps' | undefined }): Primitive {
    const { vertices, indices } = Icosahedron({ subset: options?.subset });
    const sphere = Polyhedron(vertices, indices, { detail, radius: 1 });
    if (options?.subset !== undefined) {
        sphere.normals = sphere.vertices; // prevent ugly surface of sphere subsets when merged into a full sphere
    }
    return sphere;
}
