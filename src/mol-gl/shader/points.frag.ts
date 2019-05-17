/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

export default `
precision highp float;
precision highp int;

#include common
#include common_frag_params
#include color_frag_params

#ifdef dPointFilledCircle
    uniform float uPointEdgeBleach;
#endif

const vec2 center = vec2(0.5);
const float radius = 0.5;

void main(){
    #include assign_material_color

    #if defined(dColorType_objectPicking) || defined(dColorType_instancePicking) || defined(dColorType_groupPicking)
        if (uAlpha < uPickingAlphaThreshold)
            discard; // ignore so the element below can be picked
        gl_FragColor = material;
    #elif defined(dColorType_depth)
        gl_FragColor = material;
    #else
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
`