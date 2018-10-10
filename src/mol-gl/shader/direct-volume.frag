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
uniform mat4 uModelView;
uniform mat4 uInvModelView;
uniform float uIsoValue;
uniform vec3 uBboxMin;
uniform vec3 uBboxMax;
uniform vec3 uBboxSize;
uniform vec3 uGridDim;
uniform sampler2D tTransferTex;

#if defined(dGridTexType_2d)
    uniform sampler2D tGridTex;
    uniform vec2 uGridTexDim;
#elif defined(dGridTexType_3d)
    uniform sampler3D tGridTex;
#endif

// float uIsoValue = exp(-1.5);
// float uIsoValue = 0.7;

varying vec4 vNearPos;
varying vec4 vFarPos;
varying vec3 vPosition;

#pragma glslify: transpose = require(./utils/transpose.glsl)

vec3 extractCameraPos(const in mat4 modelView) {
    // Get the 3 basis vector planes at the camera origin and transform them into model space.
    //
    // NOTE: Planes have to be transformed by the inverse transpose of a matrix
    //       Nice reference here: http://www.opengl.org/discussion_boards/showthread.php/159564-Clever-way-to-transform-plane-by-matrix
    //
    //       So for a transform to model space we need to do:
    //            inverse(transpose(inverse(MV)))
    //       This equals : transpose(MV) - see Lemma 5 in http://mathrefresher.blogspot.com.au/2007/06/transpose-of-matrix.html
    //
    // As each plane is simply (1,0,0,0), (0,1,0,0), (0,0,1,0) we can pull the data directly from the transpose matrix.
    //
    mat4 modelViewT = transpose(modelView);

    // Get plane normals
    vec3 n1 = vec3(modelViewT[0]);
    vec3 n2 = vec3(modelViewT[1]);
    vec3 n3 = vec3(modelViewT[2]);

    // Get plane distances
    float d1 = modelViewT[0].w;
    float d2 = modelViewT[1].w;
    float d3 = modelViewT[2].w;

    // Get the intersection of these 3 planes
    // (uisng math from RealTime Collision Detection by Christer Ericson)
    vec3 n2n3 = cross(n2, n3);
    float denom = dot(n1, n2n3);

    vec3 top = (n2n3 * d1) + cross(n1, (d3 * n2) - (d2 * n3));
    return top / -denom;
}

vec3 palette(in float t, in vec3 a, in vec3 b, in vec3 c, in vec3 d) {
    return a + b * cos(6.28318 * (c * t + d));
}

vec3 palette1(in float t) {
    return palette(t, vec3(0.5,0.5,0.5),vec3(0.5,0.5,0.5),vec3(1.0,1.0,1.0),vec3(0.0,0.10,0.20));
}

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

vec4 raymarch(vec3 cameraPos) {
    vec3 pos = unitCoord;
    float prevValue = -127.0;
    float value = 0.0;
    // float MAX_STEPS_F = max(max(uGridDim.x, uGridDim.y), uGridDim.z);
    // int MAX_STEPS = 2 * int(length(vec3(imgresx, imgresy, imgresz)));
    float stepSize = 1.0 / float(dMaxSteps);
    vec4 src = vec4(0.0);
    vec4 dst = vec4(0.0);

    vec3 rayDir = normalize(origPos - cameraPos);
    // rayDir = normalize(vec3(1.0, 1.0, 0.0));
    // return vec4(rayDir, 0.5);
    vec3 isoPos;
    float tmp;
    vec3 gradient = vec3(1.0);
    vec3 step = rayDir * (1.0 / uGridDim) * 0.5;

    vec3 scaleVol = vec3(1.0) / uGridDim;
    vec3 dx = vec3(gradOffset * scaleVol.x, 0.0, 0.0);
    vec3 dy = vec3(0.0, gradOffset * scaleVol.y, 0.0);
    vec3 dz = vec3(0.0, 0.0, gradOffset * scaleVol.z);

    // dst = vec4(textureVal(vec3(pos.xy, 0.6)).xyz, 0.5);
    // vec2 foo = (vec2(5.0 * uGridDim.x, 5.0 * uGridDim.y) + (pos.xy * uGridDim.xy)) / uGridTexDim;
    // dst = texture2D(tGridTex, foo);
    // dst = texture2D(tGridTex, unitCoord.xy);
    // dst.xyz = pos;
    // return mix(dst, vec4(1.0, 0.0, 0.0, 1.0), 0.5);

    for(int i = 0; i < dMaxSteps; ++i){
        if( pos.x <= 1.0 && pos.y <= 1.0 && pos.z <= 1.0 && pos.x >= 0.0 && pos.y >= 0.0 && pos.z >= 0.0) {
            value = textureVal(pos).a; // current voxel value
        } else {
            break;
        }

        #if defined(dRenderMode_volume)
            // src = texture1D(transferRGBASampler, scalarData);
            src = transferFunction(value);
            // src.rgb = palette1(value);
            // src.a = 1.0 - pow(1.0 - src.a, 0.5);
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

                float d = float(dot(gradient, normalize(cameraPos)) > 0.0);
                gradient = (2.0 * d - 1.0) * gradient;

                src.rgb = color.rgb * abs(dot(gradient, normalize(cameraPos)));
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
    // vec3 cameraPos = extractCameraPos(uModelView);
    // vec3 cameraPos = vec3(10.0, 0.0, 0.0);
    gl_FragColor = raymarch(cameraPos);
    if (length(gl_FragColor.rgb) < 0.00001) discard;
    #if defined(dRenderMode_volume)
        gl_FragColor.a = uAlpha;
    #endif
    // gl_FragColor = vec4(unitCoord, 1.0);
    // gl_FragColor = vec4(1.0, 0.0, 0.0, 0.5);
}


// const float relativeStepSize = 1.0;
// vec3 u_size = uGridDim;

// void main () {
//     // Normalize clipping plane info
//     vec3 farPos = vFarPos.xyz / vFarPos.w;
//     vec3 nearPos = vNearPos.xyz / vNearPos.w;
//     // Calculate unit vector pointing in the view direction through this fragment.
//     vec3 viewRay = normalize(nearPos.xyz - farPos.xyz);

//     // Compute the (negative) distance to the front surface or near clipping plane.
//     // v_position is the back face of the cuboid, so the initial distance calculated in the dot
//     // product below is the distance from near clip plane to the back of the cuboid
//     float distance = dot(nearPos - vPosition, viewRay);
//     distance = max(distance, min((-0.5 - vPosition.x) / viewRay.x, (u_size.x - 0.5 - vPosition.x) / viewRay.x));
//     distance = max(distance, min((-0.5 - vPosition.y) / viewRay.y, (u_size.y - 0.5 - vPosition.y) / viewRay.y));
//     distance = max(distance, min((-0.5 - vPosition.z) / viewRay.z, (u_size.z - 0.5 - vPosition.z) / viewRay.z));
//     // Now we have the starting position on the front surface
//     vec3 front = vPosition + viewRay * distance;
//     // Decide how many steps to take
//     int nsteps = int(-distance / relativeStepSize + 0.5);
//     // if (nsteps < 1)
//         // discard;
//     // Get starting location and step vector in texture coordinates
//     vec3 step = ((vPosition - front) / u_size) / float(nsteps);
//     vec3 startLoc = front / u_size;
//     // For testing: show the number of steps. This helps to establish
//     // whether the rays are correctly oriented
//     gl_FragColor = vec4(0.0, float(nsteps) / 1.0 / u_size.x, 1.0, 1.0);
//     gl_FragColor = vec4(1.0, 0.0, 0.0, 1.0);
//     return;

//     // if (u_renderstyle == 0)
//         // cast_mip(startLoc, step, nsteps, viewRay);
//     // else if (u_renderstyle == 1)
//     //     cast_iso(start_loc, step, nsteps, view_ray);


//     if (gl_FragColor.a < 0.05)
//         discard;
// }