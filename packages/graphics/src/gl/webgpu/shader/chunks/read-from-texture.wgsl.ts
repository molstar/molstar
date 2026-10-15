/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

/** Item `i` of a data texture, stored row by row */
export const read_from_texture_wgsl = /* wgsl */ `
fn readFromTexture(t: texture_2d<f32>, i: i32) -> vec4<f32> {
    let width = i32(textureDimensions(t).x);
    return textureLoad(t, vec2<i32>(i % width, i / width), 0);
}
`;

/**
 * Index of the item of the data texture `type` (a define like `dOverpaintType`) as `name`, for
 * the per-instance types; expects `instanceIndex`, `groupInstanceIndex` and `vertexInstanceIndex`
 */
export function assign_data_index_wgsl(type: string, name: string) {
  return /* wgsl */ `
#if defined(${type}_instance)
    let ${name} = instanceIndex;
#elif defined(${type}_groupInstance)
    let ${name} = groupInstanceIndex;
#elif defined(${type}_vertexInstance)
    let ${name} = vertexInstanceIndex;
#endif
`;
}
