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
 * Index of the item of the data texture `type` (a define like `dOverpaintType`, or none: the marker, clipping) as
 * `name`, for the per-instance types: 'texture' per instance or group instance (`uInstanceGranularity`),
 * 'vertexInstance'; expects `instanceIndex`, `groupInstanceIndex` and `vertexInstanceIndex`
 */
export function assign_data_index_wgsl(type: string | undefined, name: string) {
  const granularity = `let ${name} = select(groupInstanceIndex, instanceIndex, uObject.uInstanceGranularity != 0u);`;
  if (!type) return granularity;
  return /* wgsl */ `
#if defined(${type}_texture)
    ${granularity}
#elif defined(${type}_vertexInstance)
    let ${name} = vertexInstanceIndex;
#endif
`;
}
