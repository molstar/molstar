/**
 * Copyright (c) 2018-2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
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

#ifdef dPointFilledCircle
    uniform float uPointEdgeBleach;
#endif

const vec2 center = vec2(0.5);
const float radius = 0.5;

void main(){
    #include clip_pixel
    #include assign_material_color

    #if defined(dRenderVariant_pick)
        #include check_picking_alpha
        gl_FragColor = material;
    #elif defined(dRenderVariant_depth)
        gl_FragColor = material;
    #elif defined(dRenderVariant_color)
        gl_FragColor = material;

        #ifdef dPointFilledCircle
            float dist = distance(gl_PointCoord, center);
            float alpha = 1.0 - smoothstep(radius - uPointEdgeBleach, radius, dist);
            if (alpha < 0.0001) discard;
            gl_FragColor.a *= alpha;
        #endif

        #include apply_marker_color
        #include apply_fog
    #endif
}
`;