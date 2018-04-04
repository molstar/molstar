/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

// adapted from three.js, MIT License Copyright 2010-2018 three.js authors

import { Vec3 } from 'mol-math/linear-algebra'
import { computeVertexNormals, appplyRadius } from '../util'

export default function Polyhedron(_vertices: Helpers.NumberArray, _indices: Helpers.NumberArray, radius: number, detail: number) {
    radius = radius || 1;
    detail = detail || 0;

    const vertices: number[] = [];
    const indices: number[] = [];

    // the subdivision creates the vertex buffer data
    subdivide(detail);

    // all vertices should lie on a conceptual sphere with a given radius
    appplyRadius(vertices, radius);

    const normals = new Float32Array(vertices.length);
    computeVertexNormals(vertices, normals)
    // this.normalizeNormals(); // smooth normals

    return { vertices, indices, normals }

    // helper functions

    function subdivide(detail: number) {
        const a = Vec3.zero()
        const b = Vec3.zero()
        const c = Vec3.zero()

        // iterate over all faces and apply a subdivison with the given detail value
        for (let i = 0; i < _indices.length; i += 3) {
            // get the vertices of the face
            Vec3.fromArray(a, _vertices, _indices[ i + 0 ] * 3)
            Vec3.fromArray(b, _vertices, _indices[ i + 1 ] * 3)
            Vec3.fromArray(c, _vertices, _indices[ i + 2 ] * 3)

            // perform subdivision
            subdivideFace(a, b, c, detail)
        }

    }

    function subdivideFace(a: Vec3, b: Vec3, c: Vec3, detail: number) {
        const cols = Math.pow(2, detail)

        // we use this multidimensional array as a data structure for creating the subdivision
        const v: Vec3[][] = []

        // construct all of the vertices for this subdivision
        for (let i = 0; i <= cols; ++i) {
            v[i] = []

            const aj = Vec3.zero()
            Vec3.lerp(aj, a, c, i / cols)

            const bj = Vec3.zero()
            Vec3.lerp(bj, b, c, i / cols)

            const rows = cols - i
            for (let j = 0; j <= rows; ++j) {
                if (j === 0 && i === cols) {
                    v[i][j] = aj
                } else {
                    const abj = Vec3.zero()
                    Vec3.lerp(abj, aj, bj, j / rows)

                    v[i][j] = abj
                }
            }
        }

        // construct all of the faces
        for (let i = 0; i < cols; ++i) {
            for (let j = 0; j < 2 * (cols - i) - 1; ++j) {
                const k = Math.floor(j / 2)
                if (j % 2 === 0) {
                    vertices.push(...v[i][k + 1], ...v[i + 1][k], ...v[i][k])
                } else {
                    vertices.push(...v[i][k + 1], ...v[i + 1][k + 1], ...v[i + 1][k])
                }
            }
        }

        console.log(v)
    }
}