/**
 * Copyright (c) 2018-2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

// adapted from three.js, MIT License Copyright 2010-2018 three.js authors

import { Vec3 } from '../../mol-math/linear-algebra';
import { computeIndexedVertexNormals, appplyRadius } from '../util';
import { Primitive } from './primitive';

export const DefaultPolyhedronProps = {
    radius: 1,
    detail: 0
};
export type PolyhedronProps = Partial<typeof DefaultPolyhedronProps>

export function Polyhedron(_vertices: ArrayLike<number>, _indices: ArrayLike<number>, props?: PolyhedronProps): Primitive {
    const { radius, detail } = { ...DefaultPolyhedronProps, ...props };
    const builder = createBuilder();
    const { vertices, indices } = builder;

    // the subdivision creates the vertex buffer data
    subdivide(detail);

    // all vertices should lie on a conceptual sphere with a given radius
    appplyRadius(vertices, radius);

    const normals = new Float32Array(vertices.length);
    computeIndexedVertexNormals(vertices, indices, normals, vertices.length / 3, indices.length / 3);

    return {
        vertices: new Float32Array(vertices),
        normals: new Float32Array(normals),
        indices: new Uint32Array(indices)
    };

    // helper functions

    function subdivide(detail: number) {
        const a = Vec3();
        const b = Vec3();
        const c = Vec3();

        // iterate over all faces and apply a subdivison with the given detail value
        for (let i = 0; i < _indices.length; i += 3) {
            // get the vertices of the face
            Vec3.fromArray(a, _vertices, _indices[ i + 0 ] * 3);
            Vec3.fromArray(b, _vertices, _indices[ i + 1 ] * 3);
            Vec3.fromArray(c, _vertices, _indices[ i + 2 ] * 3);

            // perform subdivision
            subdivideFace(a, b, c, detail);
        }

    }

    function subdivideFace(a: Vec3, b: Vec3, c: Vec3, detail: number) {
        const cols = Math.pow(2, detail);

        // we use this multidimensional array as a data structure for creating the subdivision
        const v: Vec3[][] = [];

        // construct all of the vertices for this subdivision
        for (let i = 0; i <= cols; ++i) {
            v[i] = [];

            const aj = Vec3();
            Vec3.lerp(aj, a, c, i / cols);

            const bj = Vec3();
            Vec3.lerp(bj, b, c, i / cols);

            const rows = cols - i;
            for (let j = 0; j <= rows; ++j) {
                if (j === 0 && i === cols) {
                    v[i][j] = aj;
                } else {
                    const abj = Vec3();
                    Vec3.lerp(abj, aj, bj, j / rows);

                    v[i][j] = abj;
                }
            }
        }

        // construct all of the faces
        for (let i = 0; i < cols; ++i) {
            for (let j = 0; j < 2 * (cols - i) - 1; ++j) {
                const k = Math.floor(j / 2);
                if (j % 2 === 0) {
                    builder.add(v[i][k + 1], v[i + 1][k], v[i][k]);
                } else {
                    builder.add(v[i][k + 1], v[i + 1][k + 1], v[i + 1][k]);
                }
            }
        }
    }
}

interface Builder {
    vertices: number[]
    indices: number[]
    add: (v1: Vec3, v2: Vec3, v3: Vec3) => void
}

function createBuilder(): Builder {
    const vertices: number[] = [];
    const indices: number[] = [];

    const vertexMap = new Map<string, number>();

    function addVertex(v: Vec3) {
        const key = `${v[0].toFixed(5)}|${v[1].toFixed(5)}|${v[2].toFixed(5)}`;
        let idx = vertexMap.get(key);
        if (idx === undefined) {
            idx = vertices.length / 3;
            vertexMap.set(key, idx);
            vertices.push(...v);
        }
        return idx;
    }

    return {
        vertices,
        indices,
        add: (v1: Vec3, v2: Vec3, v3: Vec3) => {
            indices.push(addVertex(v1), addVertex(v2), addVertex(v3));
        }
    };
}