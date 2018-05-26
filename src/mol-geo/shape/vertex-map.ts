/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { ChunkedArray } from 'mol-data/util';
import { Mesh } from './mesh';

/** Mapping between vertices and ids */
interface VertexMap {
    idCount: number,
    offsetCount: number,
    ids: Helpers.NumberArray
    offsets: Uint32Array,
}

function createOffsets(idCount: number, ids: Helpers.NumberArray) {
    const offsets = ChunkedArray.create(Uint32Array, 1, 1024, 2048);
    let prevId = ids[0]
    ChunkedArray.add(offsets, 0)
    for (let i = 1; i < idCount; ++i) {
        if (prevId !== ids[i]) {
            prevId = ids[i]
            ChunkedArray.add(offsets, i)
        }
    }
    ChunkedArray.add(offsets, idCount)
    return ChunkedArray.compact(offsets, false) as Uint32Array
}

namespace VertexMap {
    export function create(idCount: number, offsetCount: number, ids: Helpers.NumberArray, offsets: Uint32Array): VertexMap {
        return {
            idCount,
            offsetCount,
            ids,
            offsets
        }
    }

    export function fromMesh(mesh: Mesh) {
        const ids = mesh.idBuffer.ref.value
        const offsets = createOffsets(mesh.vertexCount, ids)
        return create(mesh.vertexCount, offsets.length, ids, offsets)
    }

    export function rangeFromId (id: number, vertexMap: VertexMap) {
        return [0, 0]
    }
}

export default VertexMap