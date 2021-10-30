/**
 * Copyright (c) 2018-2021 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

export const points_frag = `
precision highp float;
precision highp int;

#include common
#include common_frag_params
#include color_frag_params
#include common_clip

const vec2 center = vec2(0.5);
const float radius = 0.5;

void main(){
    #include clip_pixel

    float fragmentDepth = gl_FragCoord.z;
    bool interior = false;
    #include assign_material_color

    #if defined(dPointStyle_circle)
        float dist = distance(gl_PointCoord, center);
        if (dist > radius) discard;
    #elif defined(dPointStyle_fuzzy)
        float dist = distance(gl_PointCoord, center);
        float fuzzyAlpha = 1.0 - smoothstep(0.0, radius, dist);
        if (fuzzyAlpha < 0.0001) discard;
    #endif

    #if defined(dRenderVariant_pick)
        #include check_picking_alpha
        gl_FragColor = material;
    #elif defined(dRenderVariant_depth)
        gl_FragColor = material;
    #elif defined(dRenderVariant_marking)
        gl_FragColor = material;
    #elif defined(dRenderVariant_color)
        gl_FragColor = material;

        #if defined(dPointStyle_fuzzy)
            gl_FragColor.a *= fuzzyAlpha;
        #endif

        #include apply_marker_color
        #include apply_fog
        #include wboit_write
    #endif
}
`;