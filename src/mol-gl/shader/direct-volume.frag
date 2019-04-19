/**
 * Copyright (c) 2017-2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author Michael Krone <michael.krone@uni-tuebingen.de>
 */

precision highp float;

varying vec3 unitCoord;
varying vec3 origPos;
varying float instance;

uniform mat4 uInvView;
uniform float uIsoValue;
uniform vec3 uGridDim;
uniform sampler2D tTransferTex;

uniform int uObjectId;
uniform int uInstanceCount;
uniform int uGroupCount;

uniform vec3 uHighlightColor;
uniform vec3 uSelectColor;
uniform vec2 uMarkerTexDim;
uniform sampler2D tMarker;

uniform float uAlpha;
uniform float uPickingAlphaThreshold;
uniform int uPickable;

#if defined(dGridTexType_2d)
    precision highp sampler2D;
    uniform sampler2D tGridTex;
    uniform vec3 uGridTexDim;
#elif defined(dGridTexType_3d)
    precision highp sampler3D;
    uniform sampler3D tGridTex;
#endif

#if defined(dColorType_uniform)
    uniform vec3 uColor;
#elif defined(dColorType_instance) || defined(dColorType_group) || defined(dColorType_groupInstance)
    uniform vec2 uColorTexDim;
    uniform sampler2D tColor;
#endif

#pragma glslify: import('./chunks/common.glsl')
#pragma glslify: import('./chunks/light-frag-params.glsl')

#pragma glslify: readFromTexture = require(./utils/read-from-texture.glsl, intMod=intMod, intDiv=intDiv, foo=foo) // foo=foo is a workaround for a bug in glslify
#pragma glslify: encodeFloatRGB = require(./utils/encode-float-rgb.glsl)
#pragma glslify: decodeFloatRGB = require(./utils/decode-float-rgb.glsl)
#pragma glslify: texture3dFrom2dNearest = require(./utils/texture3d-from-2d-nearest.glsl, intMod=intMod, intDiv=intDiv, foo=foo) // foo=foo is a workaround for a bug in glslify
#pragma glslify: texture3dFrom2dLinear = require(./utils/texture3d-from-2d-linear.glsl, intMod=intMod, intDiv=intDiv, foo=foo) // foo=foo is a workaround for a bug in glslify

#if defined(dGridTexType_2d)
    vec4 textureVal(vec3 pos) {
        return texture3dFrom2dLinear(tGridTex, pos, uGridDim, uGridTexDim.xy);
    }
    vec4 textureGroup(vec3 pos) {
        vec3 nearestPos = floor(pos * uGridDim + 0.5) / uGridDim + 0.5 / uGridDim;
        return texture3dFrom2dNearest(tGridTex, nearestPos, uGridDim, uGridTexDim.xy);
    }
#elif defined(dGridTexType_3d)
    vec4 textureVal(vec3 pos) {
        return texture(tGridTex, pos);
    }
    vec4 textureGroup(vec3 pos) {
        return texelFetch(tGridTex, ivec3(pos * uGridDim), 0);
    }
#endif

vec4 transferFunction(float value) {
    return texture2D(tTransferTex, vec2(value, 0.0));
}

const float gradOffset = 0.5;

vec4 raymarch(vec3 startLoc, vec3 step, vec3 viewDir) {
    vec3 scaleVol = vec3(1.0) / uGridDim;
    vec3 pos = startLoc;
    float prevValue = -1.0;
    float value = 0.0;
    vec4 src = vec4(0.0);
    vec4 dst = vec4(0.0);

    #if defined(dRenderMode_isosurface)
        vec3 isoPos;
        float tmp;

        vec3 color = vec3(0.45, 0.55, 0.8);
        vec3 gradient = vec3(1.0);
        vec3 dx = vec3(gradOffset * scaleVol.x, 0.0, 0.0);
        vec3 dy = vec3(0.0, gradOffset * scaleVol.y, 0.0);
        vec3 dz = vec3(0.0, 0.0, gradOffset * scaleVol.z);
    #endif

    for(int i = 0; i < dMaxSteps; ++i){
        value = textureVal(pos).a; // current voxel value
        if(pos.x > 1.01 || pos.y > 1.01 || pos.z > 1.01 || pos.x < -0.01 || pos.y < -0.01 || pos.z < -0.01)
            break;

        #if defined(dRenderMode_volume)
            src = transferFunction(value);
            src.rgb *= src.a;
            dst = (1.0 - dst.a) * src + dst; // standard blending
        #endif

        #if defined(dRenderMode_isosurface)
            if(prevValue > 0.0 && ( // there was a prev Value
                (prevValue < uIsoValue && value > uIsoValue) || // entering isosurface
                (prevValue > uIsoValue && value < uIsoValue) // leaving isosurface
            )) {
                tmp = ((prevValue - uIsoValue) / ((prevValue - uIsoValue) - (value - uIsoValue)));
                isoPos = mix(pos - step, pos, tmp);

                #if defined(dColorType_objectPicking) || defined(dColorType_instancePicking) || defined(dColorType_groupPicking)
                    if (uAlpha < uPickingAlphaThreshold)
                        discard; // ignore so the element below can be picked
                    if (uPickable == 0)
                        return vec4(0.0, 0.0, 0.0, 1.0); // set to empty picking id
                #endif

                #if defined(dColorType_objectPicking)
                    return vec4(encodeFloatRGB(float(uObjectId)), 1.0);
                #elif defined(dColorType_instancePicking)
                    return vec4(encodeFloatRGB(instance), 1.0);
                #elif defined(dColorType_groupPicking)
                    float group = floor(decodeFloatRGB(textureGroup(isoPos).rgb) + 0.5);
                    return vec4(encodeFloatRGB(group), 1.0);
                #else
                    // compute gradient by central differences
                    gradient.x = textureVal(isoPos - dx).a - textureVal(isoPos + dx).a;
                    gradient.y = textureVal(isoPos - dy).a - textureVal(isoPos + dy).a;
                    gradient.z = textureVal(isoPos - dz).a - textureVal(isoPos + dz).a;
                    gradient = normalize(gradient);
                    float d = float(dot(gradient, viewDir) > 0.0);
                    gradient = (2.0 * d - 1.0) * gradient;

                    float group = floor(decodeFloatRGB(textureGroup(isoPos).rgb) + 0.5);

                    #if defined(dColorType_instance)
                        color = readFromTexture(tColor, instance, uColorTexDim).rgb;
                    #elif defined(dColorType_group)
                        color = readFromTexture(tColor, group, uColorTexDim).rgb;
                    #elif defined(dColorType_groupInstance)
                        color = readFromTexture(tColor, instance * float(uGroupCount) + group, uColorTexDim).rgb;
                    #endif

                    vec3 normal = normalize(gradient);
                    vec3 vViewPosition = normalize(viewDir);
                    vec4 material = vec4(color, uAlpha);
                    #pragma glslify: import('./chunks/apply-light-color.glsl')

                    float vMarker = readFromTexture(tMarker, instance * float(uGroupCount) + group, uMarkerTexDim).a;
                    #pragma glslify: import('./chunks/apply-marker-color.glsl')

                    src.rgb = gl_FragColor.rgb;
                    src.a = gl_FragColor.a;

                    // draw interior darker
                    if( (prevValue - uIsoValue) > 0.0 ) {
                        src.rgb *= 0.5;
                    }

                    src.rgb *= src.a;
                    dst = (1.0 - dst.a) * src + dst; // standard blending
                    if(dst.a >= 1.0) {
                        break;
                    }
                #endif
            }
            prevValue = value;
        #endif

        pos += step;
    }
    return dst;
}

void main () {
    vec3 cameraPos = uInvView[3].xyz / uInvView[3].w;

    vec3 rayDir = normalize(origPos - cameraPos);
    vec3 startLoc = unitCoord;
    vec3 step = rayDir * (1.0 / uGridDim) * 0.1;

    gl_FragColor = raymarch(startLoc, step, normalize(cameraPos));
    if (length(gl_FragColor.rgb) < 0.00001) discard;
    #if defined(dRenderMode_volume)
        gl_FragColor.a *= uAlpha;
    #endif
}