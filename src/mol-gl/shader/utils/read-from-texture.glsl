/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

vec4 readFromTexture (const in sampler2D tex, const in float i, const in vec2 size) {
    float x = mod(i, size.x);
    float y = floor(i / size.x);
    vec2 uv = (vec2(x, y) + 0.5) / size;
    return texture2D(tex, uv);
}
#pragma glslify: export(readFromTexture)