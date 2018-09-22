/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

precision highp float;
precision highp int;

#pragma glslify: import('./chunks/common-frag-params.glsl')
#pragma glslify: import('./chunks/color-frag-params.glsl')

#ifdef dPointFilledCircle
    uniform float uPointEdgeBleach;
#endif

const vec2 center = vec2(0.5);
const float radius = 0.5;

void main(){
    #pragma glslify: import('./chunks/assign-material-color.glsl')

    #if defined(dColorType_objectPicking) || defined(dColorType_instancePicking) || defined(dColorType_groupPicking)
        gl_FragColor = material;
    #else
        gl_FragColor = material;

        #ifdef dPointFilledCircle
            float dist = distance(gl_PointCoord, center);
            float alpha = 1.0 - smoothstep(radius - uPointEdgeBleach, radius, dist);
            if (alpha < 0.0001) discard;
            gl_FragColor.a *= alpha;
        #endif

        #pragma glslify: import('./chunks/apply-marker-color.glsl')
        #pragma glslify: import('./chunks/apply-fog.glsl')
    #endif
}