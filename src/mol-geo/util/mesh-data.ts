/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { ValueCell } from 'mol-util/value-cell'
import { Mesh } from '../mesh/mesh';

type MeshData = {
    aPosition: ValueCell<Float32Array>,
    aNormal: ValueCell<Float32Array>,
    aGroup: ValueCell<Float32Array>,
}

export function getMeshData(mesh: Mesh): MeshData {
    return {
        aPosition: mesh.vertexBuffer,
        aNormal: mesh.normalBuffer,
        aGroup: mesh.groupBuffer,
    }
}