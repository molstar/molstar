/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

#ifdef FLAT_SHADED
    #extension GL_OES_standard_derivatives : enable
#endif

precision highp float;

// uniform vec3 light_position;
uniform vec3 light_color;
uniform vec3 light_ambient;
uniform mat4 view;
uniform float alpha;

#ifndef FLAT_SHADED
    varying vec3 vNormal;
#endif
varying vec3 vViewPosition;

#pragma glslify: import('./chunks/color-frag-params.glsl')

#pragma glslify: attenuation = require(./utils/attenuation.glsl)
#pragma glslify: calculateSpecular = require(./utils/phong-specular.glsl)
#pragma glslify: calculateDiffuse = require(./utils/oren-nayar-diffuse.glsl)

const float specularScale = 0.65;
const float shininess = 100.0;
const float roughness = 5.0;
const float albedo = 0.95;

void main() {
    // material color
    #pragma glslify: import('./chunks/color-assign-material.glsl')

    // determine surface to light direction
    // vec4 lightPosition = view * vec4(light_position, 1.0);
    // vec3 lightVector = lightPosition.xyz - vViewPosition;
    vec3 lightVector = vViewPosition;

    vec3 L = normalize(lightVector); // light direction
    vec3 V = normalize(vViewPosition); // eye direction

    // surface normal
    #ifdef FLAT_SHADED
        vec3 fdx = dFdx(vViewPosition);
        vec3 fdy = dFdy(vViewPosition);
        vec3 N = -normalize(cross(fdx, fdy));
    #else
        vec3 N = -normalize(vNormal);
        #ifdef DOUBLE_SIDED
            N = N * (float(gl_FrontFacing) * 2.0 - 1.0);
        #endif
    #endif

    // compute our diffuse & specular terms
    float specular = calculateSpecular(L, V, N, shininess) * specularScale;
    vec3 diffuse = light_color * calculateDiffuse(L, V, N, roughness, albedo);
    vec3 ambient = light_ambient;

    // add the lighting
    vec3 finalColor = material * (diffuse + ambient) + specular;

    // gl_FragColor.rgb = N;
    // gl_FragColor.a = 1.0;
    // gl_FragColor.rgb = vec3(1.0, 0.0, 0.0);
    gl_FragColor.rgb = finalColor;
    gl_FragColor.a = alpha;
}