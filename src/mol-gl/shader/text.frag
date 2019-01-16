/**
 * Copyright (c) 2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

precision highp float;
precision highp int;

#pragma glslify: import('./chunks/common-frag-params.glsl')
#pragma glslify: import('./chunks/color-frag-params.glsl')

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
    #pragma glslify: import('./chunks/assign-material-color.glsl')

    if (vTexCoord.x > 1.0) {
        gl_FragColor = vec4(uBackgroundColor, uBackgroundOpacity);
    } else {
        // TODO nicer border

        // retrieve signed distance
        float sdf = texture2D(tFont, vTexCoord).a + uBorderWidth;

        // perform adaptive anti-aliasing of the edges
        float w = clamp(
            smoothness * (abs(dFdx(vTexCoord.x)) + abs(dFdy(vTexCoord.y))),
            0.0,
            0.5
        );
        float a = smoothstep(0.5 - w, 0.5 + w, sdf);

        // gamma correction for linear attenuation
        a = pow(a, 1.0 / gamma);

        if (a < 0.5) discard;
        material.a *= a;

        if (uBorderWidth > 0.0 && sdf < (0.5 + uBorderWidth)) {
            material.xyz = uBorderColor;
        }

        gl_FragColor = material;
    }

    #if defined(dColorType_objectPicking) || defined(dColorType_instancePicking) || defined(dColorType_groupPicking)
        if (uAlpha < uPickingAlphaThreshold)
            discard; // ignore so the element below can be picked
    #else
        #pragma glslify: import('./chunks/apply-marker-color.glsl')
        #pragma glslify: import('./chunks/apply-fog.glsl')
    #endif
}