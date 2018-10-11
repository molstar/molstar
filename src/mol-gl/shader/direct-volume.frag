/**
 * Copyright (c) 2017-2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author Michael Krone <michael.krone@uni-tuebingen.de>
 */

#if defined(dGridTexType_2d)
    precision mediump sampler2D;
#elif defined(dGridTexType_3d)
    precision mediump sampler3D;
#endif
precision highp float;

varying vec3 unitCoord;
varying vec3 origPos;

uniform float uAlpha;
uniform mat4 uInvView;
uniform float uIsoValue;
uniform vec3 uGridDim;
uniform sampler2D tTransferTex;

#if defined(dGridTexType_2d)
    uniform sampler2D tGridTex;
    uniform vec2 uGridTexDim;
#elif defined(dGridTexType_3d)
    uniform sampler3D tGridTex;
#endif

#if defined(dGridTexType_2d)
    // TODO workaround due to some kind of GPU bug
    float myMod(float a, float b) {
        return a - b * float(int(a) / int(b));
    }
    float myDiv(float a, float b) {
        return float(int(a) / int(b));
    }

    vec4 textureVal(vec3 pos) {
        float zSlice0 = floor(pos.z * uGridDim.z);
        float column0 = myMod(zSlice0 * uGridDim.x, uGridTexDim.x) / uGridDim.x;
        float row0 = floor(myDiv(zSlice0 * uGridDim.x, uGridTexDim.x));
        vec2 coord0 = (vec2(column0 * uGridDim.x, row0 * uGridDim.y) + (pos.xy * uGridDim.xy)) / uGridTexDim;
        vec4 color0 = texture2D(tGridTex, coord0);

        float zSlice1 = zSlice0 + 1.0;
        float column1 = myMod(zSlice1 * uGridDim.x, uGridTexDim.x) / uGridDim.x;
        float row1 = floor(myDiv(zSlice1 * uGridDim.x, uGridTexDim.x));
        vec2 coord1 = (vec2(column1 * uGridDim.x, row1 * uGridDim.y) + (pos.xy * uGridDim.xy)) / uGridTexDim;
        vec4 color1 = texture2D(tGridTex, coord1);

        float delta0 = abs((pos.z * uGridDim.z) - zSlice0);
        return mix(color0, color1, delta0);
    }
#elif defined(dGridTexType_3d)
    vec4 textureVal(vec3 pos) {
        return texture(tGridTex, pos);
    }
#endif

vec4 transferFunction(float value) {
    return texture2D(tTransferTex, vec2(value, 0.0));
}

const float gradOffset = 0.5;
const vec3 color = vec3(0.45, 0.55, 0.8);

vec4 raymarch(vec3 startLoc, vec3 step, vec3 viewDir) {
    vec3 scaleVol = vec3(1.0) / uGridDim;
    vec3 pos = startLoc + scaleVol * 0.5;
    float prevValue = -127.0;
    float value = 0.0;
    vec4 src = vec4(0.0);
    vec4 dst = vec4(0.0);

    #if defined(dRenderMode_isosurface)
        vec3 isoPos;
        float tmp;

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

                // compute gradient by central differences
                gradient.x = textureVal(isoPos - dx).a - textureVal(isoPos + dx).a;
                gradient.y = textureVal(isoPos - dy).a - textureVal(isoPos + dy).a;
                gradient.z = textureVal(isoPos - dz).a - textureVal(isoPos + dz).a;
                gradient = normalize(gradient);

                float d = float(dot(gradient, viewDir) > 0.0);
                gradient = (2.0 * d - 1.0) * gradient;

                src.rgb = color.rgb * abs(dot(gradient, viewDir));
                src.a = uAlpha;

                // draw interior darker
                if( (prevValue - uIsoValue) > 0.0 ) {
                    src.rgb *= 0.5;
                }

                src.rgb *= src.a;
                dst = (1.0 - dst.a) * src + dst; // standard blending
                if(dst.a >= 1.0) {
                    break;
                }
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
    vec3 step = rayDir * (1.0 / uGridDim) * 0.5;

    gl_FragColor = raymarch(startLoc, step, normalize(cameraPos));
    if (length(gl_FragColor.rgb) < 0.00001) discard;
    #if defined(dRenderMode_volume)
        gl_FragColor.a = uAlpha;
    #endif
}