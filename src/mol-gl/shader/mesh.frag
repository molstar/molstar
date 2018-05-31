/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

#ifdef dFlatShaded
    #extension GL_OES_standard_derivatives : enable
#endif

precision highp float;
precision highp int;

#pragma glslify: import('./chunks/common-frag-params.glsl')
#pragma glslify: import('./chunks/color-frag-params.glsl')

// uniform vec3 uLightPosition;
uniform vec3 uLightColor;
uniform vec3 uLightAmbient;
uniform mat4 uView;
uniform float uAlpha;

#ifndef dFlatShaded
    varying vec3 vNormal;
#endif
varying vec3 vViewPosition;

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

    #if defined(dColorType_objectPicking) || defined(dColorType_instancePicking) || defined(dColorType_elementPicking)
        // gl_FragColor = vec4(material.r, material.g, material.a, 1.0);
        gl_FragColor = material;
    #else
        // determine surface to light direction
        // vec4 viewLightPosition = view * vec4(lightPosition, 1.0);
        // vec3 lightVector = viewLightPosition.xyz - vViewPosition;
        vec3 lightVector = vViewPosition;

        vec3 L = normalize(lightVector); // light direction
        vec3 V = normalize(vViewPosition); // eye direction

        // surface normal
        #ifdef dFlatShaded
            vec3 fdx = dFdx(vViewPosition);
            vec3 fdy = dFdy(vViewPosition);
            vec3 N = -normalize(cross(fdx, fdy));
        #else
            vec3 N = -normalize(vNormal);
            #ifdef dDoubleSided
                N = N * (float(gl_FrontFacing) * 2.0 - 1.0);
            #endif
        #endif

        // compute our diffuse & specular terms
        float specular = calculateSpecular(L, V, N, shininess) * specularScale;
        vec3 diffuse = uLightColor * calculateDiffuse(L, V, N, roughness, albedo);
        vec3 ambient = uLightAmbient;

        // add the lighting
        vec3 finalColor = material.rgb * (diffuse + ambient) + specular;

        // gl_FragColor.rgb = N;
        // gl_FragColor.a = 1.0;
        // gl_FragColor.rgb = vec3(1.0, 0.0, 0.0);
        gl_FragColor.rgb = finalColor;
        gl_FragColor.a = uAlpha;
    #endif
}