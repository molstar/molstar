/**
 * Copyright (c) 2020-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

/** of the data layers `tColor` (layer 0) and `tOverpaint` (1), read by their accessors (direct volumes) */
export const texture3d_from_1d_trilinear = `
vec4 readDataLayer(const in int layer, const in float i) {
    if (layer == 0) return read_tColor(i);
    return read_tOverpaint(i);
}

vec4 texture3dFrom1dTrilinear(const in int layer, const in vec3 pos, const in vec3 gridDim, const in float offset) {
    float gdYZ = gridDim.z * gridDim.y;
    float gdZ = gridDim.z;
    vec3 p0 = floor(pos * gridDim);
    vec3 p1 = ceil(pos * gridDim);
    // not \`(pos * gridDim - p0) / (p1 - p0)\`: 0 / 0 on grid coordinates
    vec3 pd = pos * gridDim - p0;
    vec4 s000 = readDataLayer(layer, offset + p0.z + p0.y * gdZ + p0.x * gdYZ);
    vec4 s100 = readDataLayer(layer, offset + p0.z + p0.y * gdZ + p1.x * gdYZ);
    vec4 s001 = readDataLayer(layer, offset + p1.z + p0.y * gdZ + p0.x * gdYZ);
    vec4 s101 = readDataLayer(layer, offset + p1.z + p0.y * gdZ + p1.x * gdYZ);
    vec4 s010 = readDataLayer(layer, offset + p0.z + p1.y * gdZ + p0.x * gdYZ);
    vec4 s110 = readDataLayer(layer, offset + p0.z + p1.y * gdZ + p1.x * gdYZ);
    vec4 s011 = readDataLayer(layer, offset + p1.z + p1.y * gdZ + p0.x * gdYZ);
    vec4 s111 = readDataLayer(layer, offset + p1.z + p1.y * gdZ + p1.x * gdYZ);
    vec4 s00 = mix(s000, s100, pd.x);
    vec4 s01 = mix(s001, s101, pd.x);
    vec4 s10 = mix(s010, s110, pd.x);
    vec4 s11 = mix(s011, s111, pd.x);
    vec4 s0 = mix(s00, s10, pd.y);
    vec4 s1 = mix(s01, s11, pd.y);
    return mix(s0, s1, pd.z);
}
`;
