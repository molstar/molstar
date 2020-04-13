/**
 * Copyright (c) 2017-2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author Michael Krone <michael.krone@uni-tuebingen.de>
 */

export default `
vec4 texture3dFrom2dLinear(sampler2D tex, vec3 pos, vec3 gridDim, vec2 texDim) {
    float zSlice0 = floor(pos.z * gridDim.z);
    float column0 = intMod(zSlice0 * gridDim.x, texDim.x) / gridDim.x;
    float row0 = floor(intDiv(zSlice0 * gridDim.x, texDim.x));
    vec2 coord0 = (vec2(column0 * gridDim.x, row0 * gridDim.y) + (pos.xy * gridDim.xy)) / texDim;
    vec4 color0 = texture2D(tex, coord0);

    float zSlice1 = zSlice0 + 1.0;
    float column1 = intMod(zSlice1 * gridDim.x, texDim.x) / gridDim.x;
    float row1 = floor(intDiv(zSlice1 * gridDim.x, texDim.x));
    vec2 coord1 = (vec2(column1 * gridDim.x, row1 * gridDim.y) + (pos.xy * gridDim.xy)) / texDim;
    vec4 color1 = texture2D(tex, coord1);

    float delta0 = abs((pos.z * gridDim.z) - zSlice0);
    return mix(color0, color1, delta0);
}
`;