/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { ValueCell } from 'mol-util/value-cell'
import { Mesh } from '../shape/mesh';

type MeshData = {
    aPosition: ValueCell<Float32Array>,
    aNormal: ValueCell<Float32Array>,
    aElementId: ValueCell<Float32Array>,
}

export function getMeshData(mesh: Mesh): MeshData {
    return {
        aPosition: mesh.vertexBuffer,
        aNormal: mesh.normalBuffer,
        aElementId: mesh.idBuffer,
    }
}