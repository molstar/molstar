/**
 * Copyright (c) 2020-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { common } from './chunks/common.glsl.js';
import { read_from_texture } from './chunks/read-from-texture.glsl.js';
import { common_vert_params } from './chunks/common-vert-params.glsl.js';
import { color_vert_params } from './chunks/color-vert-params.glsl.js';
import { common_clip } from './chunks/common-clip.glsl.js';
import { common_animation } from './chunks/common-animation.glsl.js';
import { assign_color_varying } from './chunks/assign-color-varying.glsl.js';
import { assign_marker_varying } from './chunks/assign-marker-varying.glsl.js';
import { assign_clipping_varying } from './chunks/assign-clipping-varying.glsl.js';
import { assign_size } from './chunks/assign-size.glsl.js';
import { clip_instance } from './chunks/clip-instance.glsl.js';

export const cylinders_vert = `
precision highp float;
precision highp int;

${common}
${read_from_texture}
${common_vert_params}
${color_vert_params}
${common_clip}
${common_animation}

varying mat4 vTransform;
// the cylinder's start and size, end and cap, the same at every corner of its box: two flat varyings
#if __VERSION__ == 100
    varying vec4 vCylinderStart;
    varying vec4 vCylinderEnd;
#else
    flat out vec4 vCylinderStart;
    flat out vec4 vCylinderEnd;
#endif
#if defined(dClipVariant_pixel) && defined(dClipObjects)
    // 1 where the clip objects may cut the cylinder (\`clipSphereClass\` of its bounding sphere): the clip tests skipped if 0
    #if __VERSION__ == 100
        varying float vClipSome;
    #else
        flat out float vClipSome;
    #endif
#endif

void main() {
    #if defined(dInstanceData)
        readInstanceData();
    #endif

    // the cylinder of the vertex (3 texels of \`tCylinder\`) and the corner of its box (6 vertices per cylinder)
    int item = VertexID / 6;
    #if __VERSION__ == 100
        int corner = imod(VertexID, 6);
    #else
        int corner = VertexID % 6;
    #endif
    // the data vertex ('vertex' color, size): the start for the corners at the start, else the end
    int vertexId = item * 2 + (corner < 2 ? 0 : 1);
    // start xyz and group, end xyz and radius scale, cap and color mode
    vec4 cylinderStart, cylinderEnd, cylinderFlags;
    read_tCylinder3(item * 3, cylinderStart, cylinderEnd, cylinderFlags);
    float group = cylinderStart.w;
    float colorMode = cylinderFlags.y;
    // -1 or 1 along the axis, across it and up
    vec3 mapping = vec3(
        corner < 2 ? -1.0 : 1.0,
        corner == 0 || corner == 2 || corner == 3 ? 1.0 : -1.0,
        corner == 3 || corner == 5 ? 1.0 : -1.0
    );

    ${assign_color_varying}
    ${assign_marker_varying}
    ${assign_clipping_varying}
    ${assign_size}

    mat4 transform = uObjectTransform * applyTumble(aTransform);
    vec3 wigStart = applyWiggle(cylinderStart.xyz, group, aInstance);
    vec3 wigEnd = applyWiggle(cylinderEnd.xyz, group, aInstance);
    mat4 modelTransform = uModel * transform;

    vTransform = modelTransform;
    vCylinderStart.xyz = (modelTransform * vec4(wigStart, 1.0)).xyz;
    vCylinderEnd.xyz = (modelTransform * vec4(wigEnd, 1.0)).xyz;
    vCylinderStart.w = size * cylinderEnd.w * uModelScale;
    vCylinderEnd.w = cylinderFlags.x;

    vModelPosition = (vCylinderStart.xyz + vCylinderEnd.xyz) * 0.5;

    #if defined(dClipVariant_pixel) && defined(dClipObjects)
        // pixel clipping by the bounding sphere: a cylinder clipped whole moved out of [ -w, +w ], one not cut not
        // tested per fragment
        float halfLength = length(vCylinderEnd.xyz - vCylinderStart.xyz) * 0.5;
        float boundingRadius = sqrt(halfLength * halfLength + vCylinderStart.w * vCylinderStart.w);
        // the ghost pass (\`uClipGhost\`): a cylinder whole by its center, kept where clipped, not tested per fragment
        int clipClass = uClipGhost
            ? (clipTest(vModelPosition / uModelScale) ? 0 : 2)
            : clipSphereClass(vModelPosition / uModelScale, boundingRadius / uModelScale);
        vClipSome = clipClass == 1 ? 1.0 : 0.0;
    #endif
    vec3 camDir = -mix(normalize(vModelPosition - uCameraPosition), uCameraDir, uIsOrtho);
    vec3 dir = vCylinderEnd.xyz - vCylinderStart.xyz;
    float f = mapping.x > 0.0 ? 1.0 : 0.0;
    // ensure cylinder 'dir' is pointing towards the camera
    if(dot(camDir, dir) < 0.0) {
        dir = -dir;
        f = 1.0 - f;
    }

    vec3 left = cross(camDir, dir);
    vec3 up = cross(left, dir);
    left = vCylinderStart.w * normalize(left);
    up = vCylinderStart.w * normalize(up);

    // move vertex in object-space from center to corner
    vModelPosition += mapping.x * dir + mapping.y * left + mapping.z * up;

    vec4 mvPosition = uView * vec4(vModelPosition, 1.0);
    vViewPosition = mvPosition.xyz;
    gl_Position = uProjection * mvPosition;

    // in front of the near plane (clip z is reversed): onto it, to avoid clipping the box; the fragment depth is the
    // impostor's
    if (gl_Position.z > gl_Position.w) {
        gl_Position.z = gl_Position.w;
    }

    #if defined(dDualColor) && (defined(dRenderVariant_color) || defined(dRenderVariant_tracing) || defined(dRenderVariant_emissive)) && defined(dColorType_texture)
        // dual-color mixing (group color types)
        // - for colorMode between 0 and 1 use colorMode to interpolate
        // - for colorMode == 2 do nothing, i.e., use vColor
        // - for colorMode == 3 use position on cylinder axis to interpolate
        if (dualColor) {
            if (colorMode <= 1.0) {
                vColor.rgb = mix(vColor.rgb, color2.rgb, colorMode);
            } else if (colorMode == 3.0) {
                vColor.rgb = mix(vColor.rgb, color2.rgb, mix(-0.25, 1.25, f / 1.5));
            }
        }
    #endif

    #if defined(dClipVariant_pixel) && defined(dClipObjects)
        if (clipClass == 2) {
            gl_Position.z = 2.0 * gl_Position.w;
        }
    #endif

    ${clip_instance}
}
`;
