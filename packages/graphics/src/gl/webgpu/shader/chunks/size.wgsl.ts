/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

/**
 * Sizes from the size data, free of any binding: shared by the render stages (`assign_size_wgsl`) and compute kernels
 * that size the same items; expects `unpackRGBToInt` (`pack_int_wgsl`)
 */
export const size_wgsl = /* wgsl */ `
/**
 * The index into the size data of an item by the size type's \`indexType\` (\`SizeTypeIndex\`): instance, group,
 * groupInstance, vertex, vertexInstance
 */
fn sizeIndex(indexType: i32, instance: i32, group: i32, groupInstance: i32, vertex: i32, vertexInstance: i32) -> i32 {
    switch indexType {
        case 1: { return group; }
        case 2: { return groupInstance; }
        case 3: { return vertex; }
        case 4: { return vertexInstance; }
        default: { return instance; }
    }
}

/** A size from its element of the size data (\`read_tSize\`): stored times 100 */
fn unpackSize(packed: vec4<f32>) -> f32 {
    return unpackRGBToInt(packed.rgb) / 100.0;
}
`;
