/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { PagedCopyShift } from '../../culled-data.js';

/**
 * A pick's texel resolved (`paged-pick.ts`): of the page store's object `uObjectId` its instance word its composed
 * pose's slot, its record's store index and copy written (`outInstance`, `outCopy`); of another object (the first
 * dispatch of a pick alone, `uFirst`) its instance as picked, no copy. Ids as the pick packs them (`packPickId`: the id
 * plus one, its bytes r, g, b its low 24 bits, a its top 8)
 */
export const paged_pick_wgsl = /* wgsl */ `
fn unpackId(c: vec4<f32>) -> u32 {
    let b = vec4<u32>(round(c * 255.0));
    return ((b.a << 24u) | (b.r << 16u) | (b.g << 8u) | b.b) - 1u;
}

fn packId(id: u32) -> vec4<f32> {
    let v = id + 1u;
    return vec4<f32>(f32((v >> 16u) & 255u), f32((v >> 8u) & 255u), f32(v & 255u), f32(v >> 24u)) / 255.0;
}

@compute @workgroup_size(WorkgroupSizeX, WorkgroupSizeY, WorkgroupSizeZ)
fn main(@builtin(global_invocation_id) id: vec3<u32>) {
    let size = vec2<u32>(uParams.uSize);
    if (any(id.xy >= size)) {
        return;
    }
    let xy = vec2<i32>(id.xy) + uParams.uOffset;
    let object = textureLoad(tObject, xy, 0);
    let instance = textureLoad(tInstance, xy, 0);
    // nothing picked (packed 0): its id 2^32 - 1
    let objectId = unpackId(object);
    if (objectId != u32(uParams.uObjectId)) {
        if (uParams.uFirst != 0) {
            textureStore(outInstance, xy, instance);
            textureStore(outCopy, xy, vec4<f32>(0.0));
        }
        return;
    }
    // the slot's record: its second vec4 (the quaternion, the store index, the segment and the copy plus one)
    let slot = unpackId(instance);
    let words = poses[poses[0].z + 2u * slot + 1u];
    textureStore(outInstance, xy, packId(words.z));
    textureStore(outCopy, xy, packId(words.w >> ${PagedCopyShift}u));
}
`;
