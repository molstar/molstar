/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Mesh } from './mesh';

/** Mapping between vertices and ids */
interface VertexMap {
    idCount: number,
    offsetCount: number,
    ids: Helpers.NumberArray | undefined
    offsets: Helpers.NumberArray,
}

function createOffsets(ids: Helpers.NumberArray | undefined) {
    return []
}

namespace VertexMap {
    export function create(idCount: number, offsetCount: number, ids: Helpers.NumberArray | undefined, offsets: Helpers.NumberArray): VertexMap {
        return {
            idCount,
            offsetCount,
            ids,
            offsets
        }
    }

    export function fromMesh(mesh: Mesh) {
        const ids = mesh.idBuffer.ref.value
        const offsets = createOffsets(ids)
        return create(mesh.vertexCount, offsets.length, ids, offsets)
    }

    export function rangeFromId (id: number, vertexMap: VertexMap) {
        return [0, 0]
    }
}

export default VertexMap