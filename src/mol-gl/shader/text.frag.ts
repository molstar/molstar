/**
 * Copyright (c) 2019-2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

export default `
precision highp float;
precision highp int;

#include common
#include common_frag_params
#include color_frag_params
#include common_clip

uniform sampler2D tFont;

uniform vec3 uBorderColor;
uniform float uBorderWidth;
uniform vec3 uBackgroundColor;
uniform float uBackgroundOpacity;

varying vec2 vTexCoord;

const float smoothness = 32.0;
const float gamma = 2.2;

void main2(){
    gl_FragColor = vec4(1.0, 0.0, 0.0, 1.0);
}

void main(){
    #include clip_pixel
    #include assign_material_color

    if (vTexCoord.x > 1.0) {
        gl_FragColor = vec4(uBackgroundColor, uBackgroundOpacity);
    } else {
        // retrieve signed distance
        float sdf = texture2D(tFont, vTexCoord).a + uBorderWidth;

        // perform adaptive anti-aliasing of the edges
        float w = clamp(smoothness * (abs(dFdx(vTexCoord.x)) + abs(dFdy(vTexCoord.y))), 0.0, 0.5);
        float a = smoothstep(0.5 - w, 0.5 + w, sdf);

        // gamma correction for linear attenuation
        a = pow(a, 1.0 / gamma);

        if (a < 0.5) discard;
        material.a *= a;

        // add border
        float t = 0.5 + uBorderWidth;
        if (uBorderWidth > 0.0 && sdf < t) {
            material.xyz = mix(uBorderColor, material.xyz, smoothstep(t - w, t, sdf));
        }

        gl_FragColor = material;
    }

    #if defined(dRenderVariant_pick)
        #include check_picking_alpha
    #elif defined(dRenderVariant_depth)
        gl_FragColor = material;
    #elif defined(dRenderVariant_color)
        #include apply_marker_color
        #include apply_fog
    #endif
}
`;