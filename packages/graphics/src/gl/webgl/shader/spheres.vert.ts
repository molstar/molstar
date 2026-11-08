/**
 * Copyright (c) 2019-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { common } from './chunks/common.glsl.js';
import { assign_quad_position_group } from './chunks/assign-quad-position-group.glsl.js';
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
import { assign_lod_factor } from './chunks/assign-lod-factor.glsl.js';

export const spheres_vert = `
precision highp float;
precision highp int;

${common}
${read_from_texture}
${common_vert_params}
${color_vert_params}
${common_clip}
${common_animation}

// the sphere's view-space center (negated) and radius, the same at every corner: one flat varying
#if __VERSION__ == 100
    varying vec4 vSphere;
#else
    flat out vec4 vSphere;
#endif
varying vec3 vPoint;

#if !defined(dClipPrimitive) && defined(dClipVariant_pixel) && defined(dClipObjects)
    #define dClipSphere
#endif
#if defined(dClipSphere)
    // 1 where the clip objects may cut the sphere (\`clipSphereClass\`): the fragment stage's clip tests skipped if 0
    #if __VERSION__ == 100
        varying float vClipSome;
    #else
        flat out float vClipSome;
    #endif
#endif

/**
 * Bounding rectangle of a clipped, perspective-projected 3D Sphere.
 * Michael Mara, Morgan McGuire. 2013
 *
 * Specialization by Arseny Kapoulkine, MIT License Copyright (c) 2018
 * https://github.com/zeux/niagara
 *
 * With the off-axis terms of asymmetric projections (stereo, XR), \`uProjection[2].xy\`, zero for symmetric ones.
 */
void sphereProjection(const in vec3 p, const in float r, const in vec2 mapping) {
    vec3 pr = p * r;
    float pzr2 = p.z * p.z - r * r;

    float vx = sqrt(p.x * p.x + pzr2);
    float minx = ((vx * p.x - pr.z) / (vx * p.z + pr.x)) * uProjection[0][0];
    float maxx = ((vx * p.x + pr.z) / (vx * p.z - pr.x)) * uProjection[0][0];

    float vy = sqrt(p.y * p.y + pzr2);
    float miny = ((vy * p.y - pr.z) / (vy * p.z + pr.y)) * uProjection[1][1];
    float maxy = ((vy * p.y + pr.z) / (vy * p.z - pr.y)) * uProjection[1][1];

    gl_Position.xy = vec2(maxx + minx, maxy + miny) * -0.5;
    gl_Position.xy -= mapping * vec2(maxx - minx, maxy - miny) * 0.5;
    gl_Position.xy -= uProjection[2].xy;
    gl_Position.xy *= gl_Position.w;
}

// the least and greatest x / z over the part beyond the near plane (z <= -near) of the circle at (a, c) of radius r (a
// sphere's projection onto the plane of an axis and z): at the tangent points from the eye that lie beyond it, and at
// the ends of the chord the plane cuts (the only bounds where the eye is in the circle); least above greatest where no
// part lies beyond it. After "2D Polyhedral Bounds of a Clipped, Perspective-Projected 3D Sphere", Michael Mara, Morgan
// McGuire. 2013
vec2 nearAxisBounds(const in float a, const in float c, const in float r, const in float near) {
    float lo = 1e30;
    float hi = -1e30;
    float l2 = a * a + c * c;
    float t2 = l2 - r * r;
    if (t2 > 0.0) {
        float t = sqrt(t2);
        for (int i = 0; i < 2; ++i) {
            float s = i == 0 ? -1.0 : 1.0;
            float tx = (a * t2 - s * c * t * r) / l2;
            float tz = (c * t2 + s * a * t * r) / l2;
            if (tz <= -near) {
                lo = min(lo, tx / tz);
                hi = max(hi, tx / tz);
            }
        }
    }
    float d = c + near;
    float w2 = r * r - d * d;
    if (w2 > 0.0) {
        float w = sqrt(w2);
        lo = min(lo, min((a - w) / -near, (a + w) / -near));
        hi = max(hi, max((a - w) / -near, (a + w) / -near));
    }
    return vec2(lo, hi);
}

// the quad of a sphere through the near plane: the bounds of its part beyond it (the whole sphere's are unbounded where
// it reaches behind the eye, a quad over the view), at a depth within the clip volume (the fragments write theirs; the
// corners on one view plane, so the rays between them interpolate exactly); false where none of it lies beyond
bool sphereNearProjection(const in vec3 p, const in float r, const in vec2 mapping) {
    vec2 x = nearAxisBounds(p.x, p.z, r, uNear) * uProjection[0][0];
    vec2 y = nearAxisBounds(p.y, p.z, r, uNear) * uProjection[1][1];
    gl_Position.xy = vec2(x.y + x.x, y.y + y.x) * -0.5;
    gl_Position.xy -= mapping * vec2(x.y - x.x, y.y - y.x) * 0.5;
    gl_Position.xy -= uProjection[2].xy;
    gl_Position.zw = vec2(0.5, 1.0);
    return x.x <= x.y;
}

void main(void){
    ${assign_quad_position_group}

    position = applyWiggle(position, group, aInstance);
    mat4 transform = uObjectTransform * applyTumble(aTransform);

    ${assign_color_varying}
    ${assign_marker_varying}
    ${assign_clipping_varying}
    ${assign_size}

    vSphere.w = size * uModelScale;

    vec4 position4 = vec4(position, 1.0);
    vModelPosition = (uModel * transform * position4).xyz; // for clipping in frag shader

    ${assign_lod_factor}
    vSphere.w *= lodFactor;

    vec4 mvPosition = uModelView * transform * position4;

    // through the near plane with none of it beyond
    bool nearNone = false;
    #ifdef dApproximate
        vec4 mvCorner = vec4(mvPosition.xyz, 1.0);
        mvCorner.xy += mapping * vSphere.w;
        gl_Position = uProjection * mvCorner;
    #else
        if (uIsOrtho == 1.0) {
            vec4 mvCorner = vec4(mvPosition.xyz, 1.0);
            mvCorner.xy += mapping * vSphere.w;
            gl_Position = uProjection * mvCorner;
        } else if (-mvPosition.z - vSphere.w > uNear) {
            gl_Position = uProjection * vec4(mvPosition.xyz, 1.0);
            sphereProjection(mvPosition.xyz, vSphere.w, mapping);
        } else {
            nearNone = !sphereNearProjection(mvPosition.xyz, vSphere.w, mapping);
        }
    #endif

    vec4 vPoint4 = uInvProjection * gl_Position;
    vPoint = vPoint4.xyz / vPoint4.w;
    vSphere.xyz = -mvPosition.xyz / mvPosition.w;

    // in front of the near plane (clip z is reversed)
    if (gl_Position.z > gl_Position.w) {
        mvPosition.z -= 2.0 * vSphere.w; // avoid clipping
        gl_Position.z = (uProjection * vec4(mvPosition.xyz, 1.0)).z;
    }

    if (lodCulled || nearNone) {
        // move out of [ -w, +w ] to 'discard' in vert shader
        gl_Position.z = 2.0 * gl_Position.w;
    }

    #if defined(dRenderVariant_marking)
        // the marking depth draws unmarked spheres only (\`uMarkingType\` 1), the mask marked ones: the others moved out of
        // [ -w, +w ], not shaded to be discarded
        float vertexMarker = uMarker == -1.0 ? floor(vMarker * 255.0 + 0.5) : uMarker;
        if ((uMarkingType == 1) == (vertexMarker > 0.0)) {
            gl_Position.z = 2.0 * gl_Position.w;
        }
    #endif

    #if defined(dClipPrimitive) && !defined(dClipVariant_instance) && defined(dClipObjects)
        // the ghost pass (\`uClipGhost\`) keeps what the clip objects remove
        if (clipTest(vModelPosition / uModelScale) != uClipGhost) {
            // move out of [ -w, +w ] to 'discard' in vert shader
            gl_Position.z = 2.0 * gl_Position.w;
        }
    #elif defined(dClipSphere)
        if (uClipGhost) {
            // the ghost pass: a sphere whole by its center, not tested per fragment
            if (!clipTest(vModelPosition / uModelScale)) gl_Position.z = 2.0 * gl_Position.w;
            vClipSome = 0.0;
        } else {
            // pixel clipping: a sphere clipped whole moved out of [ -w, +w ], one not cut not tested per fragment
            int clipClass = clipSphereClass(vModelPosition / uModelScale, vSphere.w / uModelScale);
            if (clipClass == 2) {
                gl_Position.z = 2.0 * gl_Position.w;
            }
            vClipSome = clipClass == 1 ? 1.0 : 0.0;
        }
    #else
        ${clip_instance}
    #endif
}
`;
