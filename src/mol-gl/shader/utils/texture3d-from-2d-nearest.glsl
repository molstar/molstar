/**
 * Copyright (c) 2017-2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author Michael Krone <michael.krone@uni-tuebingen.de>
 */

// TODO workaround due to some kind of GPU quirk
float myMod(float a, float b) { return a - b * float(int(a) / int(b)); }
float myDiv(float a, float b) { return float(int(a) / int(b)); }

vec4 texture3dFrom2dNearest(sampler2D tex, vec3 pos, vec3 gridDim, vec2 texDim) {
    float zSlice = floor(pos.z * gridDim.z + 0.5); // round to nearest z-slice
    float column = myMod(zSlice * gridDim.x, texDim.x) / gridDim.x;
    float row = floor(myDiv(zSlice * gridDim.x, texDim.x));
    vec2 coord = (vec2(column * gridDim.x, row * gridDim.y) + (pos.xy * gridDim.xy)) / texDim;
    return texture2D(tex, coord);
}

#pragma glslify: export(texture3dFrom2dNearest)