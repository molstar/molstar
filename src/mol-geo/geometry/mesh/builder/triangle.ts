/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import {  Mat4 } from 'mol-math/linear-algebra';
import { MeshBuilder } from '../mesh-builder';

const tmpSphereMat = Mat4.identity()

function getTriangle(vertices: number[], normals: number[], indices: number[]) {

    return {vertices, normals, indices};
}

export function addTriangle(state: MeshBuilder.State, triangle_vertices: number[], triangle_normals: number[], triangle_indices: number[]) {
    MeshBuilder.addPrimitive(state, tmpSphereMat, getTriangle( triangle_vertices, triangle_normals, triangle_indices))
}