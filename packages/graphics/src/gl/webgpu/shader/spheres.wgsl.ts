/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { animation_wgsl, assign_wiggle_wgsl } from './chunks/animation.wgsl.js';
import { clip_flag_wgsl, clip_functions_wgsl, clip_instance_wgsl, clip_pixel_wgsl } from './chunks/clip.wgsl.js';
import { segment_wgsl } from './chunks/segment.wgsl.js';
import { common_wgsl } from './chunks/common.wgsl.js';
import { light_wgsl } from './chunks/light.wgsl.js';
import {
  apply_interior_wgsl,
  apply_light_wgsl,
  apply_fog_wgsl,
  apply_marker_wgsl,
  assign_material_wgsl,
  check_transparency_wgsl,
  material_functions_wgsl,
} from './chunks/material.wgsl.js';
import {
  depth_output_wgsl,
  dpoit_output_wgsl,
  emissive_output_wgsl,
  fragment_output_wgsl,
  marking_output_wgsl,
  pick_output_wgsl,
  tracing_output_wgsl,
  wboit_output_wgsl,
} from './chunks/output.wgsl.js';
import { packing_wgsl } from './chunks/packing.wgsl.js';
import { read_from_texture_wgsl } from './chunks/read-from-texture.wgsl.js';
import {
  assign_quad_position_group_wgsl,
  assign_size_wgsl,
  assign_varyings_wgsl,
  vertex_output_wgsl,
} from './chunks/varyings.wgsl.js';

/** Sphere impostors: a screen-aligned quad per sphere covering its projection, ray-sphere intersection per fragment */
export const spheres_wgsl = /* wgsl */ `
${common_wgsl}
${read_from_texture_wgsl}
${animation_wgsl}
${light_wgsl}
${material_functions_wgsl}
${segment_wgsl}
${clip_functions_wgsl}
${packing_wgsl}

${vertex_output_wgsl(`    @location(9) vRadius: f32,
    @location(10) vPoint: vec3<f32>,
    @location(11) vPointViewPosition: vec3<f32>,`)}

/**
 * Bounding rectangle of a clipped, perspective-projected 3D Sphere.
 * Michael Mara, Morgan McGuire. 2013
 *
 * Specialization by Arseny Kapoulkine, MIT License Copyright (c) 2018
 * https://github.com/zeux/niagara
 *
 * Only works for for symmetric projections.
 */
fn sphereProjection(p: vec3<f32>, r: f32, mapping: vec2<f32>, w: f32) -> vec2<f32> {
    let pr = p * r;
    let pzr2 = p.z * p.z - r * r;
    let projection = uGlobal.uProjection;

    let vx = sqrt(p.x * p.x + pzr2);
    let minx = ((vx * p.x - pr.z) / (vx * p.z + pr.x)) * projection[0][0];
    let maxx = ((vx * p.x + pr.z) / (vx * p.z - pr.x)) * projection[0][0];

    let vy = sqrt(p.y * p.y + pzr2);
    let miny = ((vy * p.y - pr.z) / (vy * p.z + pr.y)) * projection[1][1];
    let maxy = ((vy * p.y + pr.z) / (vy * p.z - pr.y)) * projection[1][1];

    let xy = vec2<f32>(maxx + minx, maxy + miny) * -0.5 - mapping * vec2<f32>(maxx - minx, maxy - miny) * 0.5;
    return xy * w;
}

/**
 * Point size and center using the technique described in:
 * "GPU-Based Ray-Casting of Quadratic Surfaces" http://dl.acm.org/citation.cfm?id=2386396
 * by Christian Sigg, Tim Weyrich, Mario Botsch, Markus Gross.
 */
fn quadraticProjection(position: vec3<f32>, radius: f32, mapping: vec2<f32>, transform: mat4x4<f32>, w: f32) -> vec2<f32> {
    let T = mat4x4<f32>(
        radius, 0.0, 0.0, 0.0,
        0.0, radius, 0.0, 0.0,
        0.0, 0.0, radius, 0.0,
        position.x, position.y, position.z, 1.0
    );
    let D = mat4x4<f32>(
        1.0, 0.0, 0.0, 0.0,
        0.0, 1.0, 0.0, 0.0,
        0.0, 0.0, 1.0, 0.0,
        0.0, 0.0, 0.0, -1.0
    );
    // only the x, y and w rows of the projection are used, they are the same for [0, 1] depth
    let R = transpose(uGlobal.uProjection * uGlobal.uModelView * transform * T);

    var A = dot(R[3], D * R[3]);
    var B = -2.0 * dot(R[0], D * R[3]);
    var C = dot(R[0], D * R[0]);
    let xbc = vec2<f32>(
        (-B - sqrt(B * B - 4.0 * A * C)) / (2.0 * A),
        (-B + sqrt(B * B - 4.0 * A * C)) / (2.0 * A)
    );
    let sx = abs(xbc.x - xbc.y) * 0.5;

    A = dot(R[3], D * R[3]);
    B = -2.0 * dot(R[1], D * R[3]);
    C = dot(R[1], D * R[1]);
    let ybc = vec2<f32>(
        (-B - sqrt(B * B - 4.0 * A * C)) / (2.0 * A),
        (-B + sqrt(B * B - 4.0 * A * C)) / (2.0 * A)
    );
    let sy = abs(ybc.x - ybc.y) * 0.5;

    let xy = vec2<f32>(0.5 * (xbc.x + xbc.y), 0.5 * (ybc.x + ybc.y)) - mapping * vec2<f32>(sx, sy);
    return xy * w;
}

@vertex
fn vertexMain(in: VertexInput) -> VertexOutput {
    var out: VertexOutput;

    ${assign_quad_position_group_wgsl}
    ${assign_wiggle_wgsl}
    position = applyWiggle(position, group, wiggle);
    let aTransform = mat4x4<f32>(in.aTransform0, in.aTransform1, in.aTransform2, in.aTransform3);
    let transform = applyTumble(aTransform, in.aInstance, f32(uObject.uObjectId), invariantBoundingSphere(in));

    ${assign_varyings_wgsl}
    ${assign_size_wgsl}

    let modelScale = uGlobal.uModelScale;
    var radius = size * modelScale;

    let position4 = vec4<f32>(position, 1.0);
    out.vModelPosition = (uGlobal.uModel * transform * position4).xyz; // for clipping in the fragment stage

    let lod = uDraw.uLod;
    let hasLod = lod.w != 0.0 && (lod.x != 0.0 || lod.y != 0.0);
    var d = 0.0;
    if (hasLod) {
        #if defined(dSegmented)
            // members of a merged item can have different scales for the same lod level (sphere-count dependent)
            let lodScale = in.aSegmentLod[clamp(uDraw.uLodLevel, 0, 3)];
        #else
            let lodScale = lod.w;
        #endif
        if (modelScale != 1.0) {
            radius *= lodScale;
        } else {
            d = (dot(uGlobal.uCameraPlane.xyz, out.vModelPosition) + uGlobal.uCameraPlane.w) / modelScale;
            let f = min(smoothstep(lod.x, lod.x + lod.z, d), 1.0 - smoothstep(lod.y - lod.z, lod.y, d)) * lodScale;
            radius *= f;
        }
    }

    var mvPosition = uGlobal.uModelView * transform * position4;
    out.vViewPosition = mvPosition.xyz;

    #if defined(dApproximate)
        out.position = uGlobal.uProjection * vec4<f32>(mvPosition.xy + mapping * radius, mvPosition.z, 1.0);
    #else
        if (uGlobal.uIsOrtho == 1.0) {
            out.position = uGlobal.uProjection * vec4<f32>(mvPosition.xy + mapping * radius, mvPosition.z, 1.0);
        } else if (uGlobal.uIsAsymmetricProjection != 0u) {
            let p = uGlobal.uProjection * vec4<f32>(mvPosition.xyz, 1.0);
            out.position = vec4<f32>(quadraticProjection(position, radius / modelScale, mapping, transform, p.w), p.zw);
        } else {
            let p = uGlobal.uProjection * vec4<f32>(mvPosition.xyz, 1.0);
            out.position = vec4<f32>(sphereProjection(mvPosition.xyz, radius, mapping, p.w), p.zw);
        }
    #endif

    out.vRadius = radius;
    let point4 = uGlobal.uInvProjection * out.position;
    out.vPoint = point4.xyz / point4.w;
    out.vPointViewPosition = -mvPosition.xyz / mvPosition.w;

    // in front of the near plane (clip z is reversed)
    if (out.position.z > out.position.w) {
        mvPosition.z -= 2.0 * radius; // avoid clipping
        out.position.z = (uGlobal.uProjection * vec4<f32>(mvPosition.xyz, 1.0)).z;
    }

    if (modelScale == 1.0 && hasLod && (d < lod.x || d > lod.y)) {
        // move out of [0, w] to 'discard' in the vertex stage
        out.position.z = 2.0 * out.position.w;
    }

    #if defined(dClipPrimitive) && !defined(dClipVariant_instance) && defined(dClipObjects)
        if (clipTest(out.vModelPosition / modelScale, clipFlag)) {
            // move out of [0, w] to 'discard' in the vertex stage
            out.position.z = 2.0 * out.position.w;
        }
    #else
        ${clip_instance_wgsl}
    #endif
    return out;
}

struct SphereHit {
    hit: bool,
    modelPosition: vec3<f32>,
    viewPosition: vec3<f32>,
    normal: vec3<f32>,
    interior: bool,
    depth: f32,
    /** on a clip cap of a solid interior, not tested by pixel clipping */
    isCap: bool,
};

#if !defined(dClipPrimitive) && defined(dClipVariant_pixel) && defined(dClipObjects)
    #define dClipSphere
#endif

fn sphereImpostor(radius: f32, point: vec3<f32>, pointViewPosition: vec3<f32>, clipFlag: i32) -> SphereHit {
    var h: SphereHit;
    h.hit = false;
    h.isCap = false;

    let isOrtho = uGlobal.uIsOrtho;
    let cameraSpherePos = -pointViewPosition;
    let rayOrigin = mix(vec3<f32>(0.0), point, isOrtho);
    let rayDirection = mix(normalize(point), vec3<f32>(0.0, 0.0, 1.0), isOrtho);
    let cameraSphereDir = mix(cameraSpherePos, rayOrigin - cameraSpherePos, isOrtho);

    let B = dot(rayDirection, cameraSphereDir);
    let det = B * B + radius * radius - dot(cameraSphereDir, cameraSphereDir);
    if (det < 0.0) {
        return h;
    }

    let sqrtDet = sqrt(det);
    let posT = mix(B + sqrtDet, B - sqrtDet, isOrtho);
    let negT = mix(B - sqrtDet, B + sqrtDet, isOrtho);

    #if defined(dRenderVariant_depth)
        if (uGlobal.uDepthBack != 0u) {
            // back-depth pass wants the farthest point that is not clipped
            h.viewPosition = rayDirection * posT + rayOrigin;
            h.modelPosition = (uGlobal.uInvView * vec4<f32>(h.viewPosition, 1.0)).xyz;
            h.normal = -normalize(h.viewPosition - cameraSpherePos);
            h.interior = true;
            #if defined(dSolidInterior) && defined(dClipSphere)
                let nearModelPosition = (uGlobal.uInvView * vec4<f32>(rayDirection * negT + rayOrigin, 1.0)).xyz;
                let s = clipExit(h.modelPosition / uGlobal.uModelScale, nearModelPosition / uGlobal.uModelScale, 0.0, clipFlag);
                if (s < 0.0) {
                    return h;
                }
                if (s > 0.0) {
                    h.viewPosition = rayDirection * mix(posT, negT, s) + rayOrigin;
                    h.modelPosition = (uGlobal.uInvView * vec4<f32>(h.viewPosition, 1.0)).xyz;
                    let invView = uGlobal.uInvView;
                    let invView3 = mat3x3<f32>(invView[0].xyz, invView[1].xyz, invView[2].xyz);
                    h.normal = -normalize(clipNormal(h.modelPosition / uGlobal.uModelScale, clipFlag) * invView3);
                    h.isCap = true;
                }
            #endif
            h.depth = calcDepth(h.viewPosition);
            h.hit = h.depth < 1.0;
            return h;
        }
    #endif

    h.viewPosition = rayDirection * negT + rayOrigin;
    h.modelPosition = (uGlobal.uInvView * vec4<f32>(h.viewPosition, 1.0)).xyz;
    h.depth = calcDepth(h.viewPosition);
    let frontDepth = h.depth;
    let frontModelPosition = h.modelPosition;
    #if defined(dClipSphere)
        if (clipTest(h.modelPosition / uGlobal.uModelScale, clipFlag)) {
            h.depth = 2.0; // clipped, beyond the near plane (depth is reversed)
        }
    #endif
    if (h.depth < 1.0) {
        h.normal = normalize(h.viewPosition - cameraSpherePos);
        h.interior = false;
        h.hit = true;
        return h;
    }

    #if defined(dSolidInterior)
        let backFaces = true;
    #else
        let backFaces = uMaterial.uDoubleSided != 0u;
    #endif
    if (backFaces) {
        h.viewPosition = rayDirection * posT + rayOrigin;
        h.modelPosition = (uGlobal.uInvView * vec4<f32>(h.viewPosition, 1.0)).xyz;
        h.depth = calcDepth(h.viewPosition);
        h.normal = -normalize(h.viewPosition - cameraSpherePos);
        h.interior = true;
        h.hit = h.depth < 1.0;
        #if defined(dSolidInterior)
            if (h.hit) {
                let nearT = -(uGlobal.uNear + rayOrigin.z) / rayDirection.z;
                let sNear = select(clamp((nearT - negT) / (posT - negT), 0.0, 1.0), 0.0, frontDepth < 1.0);
                #if defined(dClipSphere)
                    // the first point along the ray inside the sphere that no clip object removes
                    let sCap = clipExit(frontModelPosition / uGlobal.uModelScale, h.modelPosition / uGlobal.uModelScale, sNear, clipFlag);
                    if (sCap < 0.0) {
                        h.hit = false;
                        return h;
                    }
                #else
                    let sCap = sNear;
                #endif
                if (sCap == sNear && frontDepth >= 1.0) {
                    // the front is clipped by the near plane: cap at the near plane
                    h.depth = 1.0 - (0.0000001 / radius);
                    h.normal = -mix(normalize(point), vec3<f32>(0.0, 0.0, -1.0), isOrtho);
                    h.viewPosition = rayDirection * nearT + rayOrigin;
                    h.modelPosition = (uGlobal.uInvView * vec4<f32>(h.viewPosition, 1.0)).xyz;
                }
                #if defined(dClipSphere)
                    else {
                        // cap on the clip object surface
                        h.viewPosition = rayDirection * mix(negT, posT, sCap) + rayOrigin;
                        h.modelPosition = (uGlobal.uInvView * vec4<f32>(h.viewPosition, 1.0)).xyz;
                        h.depth = calcDepth(h.viewPosition) - (0.0000001 / radius);
                        let invView = uGlobal.uInvView;
                        let invView3 = mat3x3<f32>(invView[0].xyz, invView[1].xyz, invView[2].xyz);
                        h.normal = -normalize(clipNormal(h.modelPosition / uGlobal.uModelScale, clipFlag) * invView3);
                        h.isCap = true;
                    }
                #endif
            }
        #endif
    }
    return h;
}

${fragment_output_wgsl(`#if !defined(dApproximate)
    @builtin(frag_depth) fragDepth: f32,
#endif`)}

@fragment
fn fragmentMain(in: VertexOutput) -> FragmentOutput {
    var out: FragmentOutput;
    ${clip_flag_wgsl('in.vClipping')}

    #if defined(dApproximate)
        var pointDir = -in.vPointViewPosition - in.vPoint;
        if (dot(pointDir, pointDir) > in.vRadius * in.vRadius) {
            discard;
        }
        let viewPosition = -in.vPointViewPosition;
        let modelPosition = in.vModelPosition;
        let fragmentDepth = in.position.z;
        pointDir.z -= cos(length(pointDir)) * in.vRadius * 0.5;
        let cameraNormal = -normalize(pointDir);
        let interior = false;
        let isCap = false;
    #else
        let h = sphereImpostor(in.vRadius, in.vPoint, in.vPointViewPosition, clipFlag);
        if (!h.hit || h.depth < 0.0 || h.depth > 1.0) {
            discard;
        }
        out.fragDepth = h.depth;
        let viewPosition = h.viewPosition;
        let modelPosition = h.modelPosition;
        let fragmentDepth = h.depth;
        let cameraNormal = h.normal;
        let interior = h.interior;
        let isCap = h.isCap;
    #endif

    if (fadeLod(in.position, modelPosition)) {
        discard;
    }
    #if !defined(dClipPrimitive)
        if (!isCap) {
            ${clip_pixel_wgsl('modelPosition')}
        }
    #endif

    var normal = -cameraNormal;

    ${assign_material_wgsl}

    // thin spheres are more transparent; marking draws all masks at once, so it checks the material alpha
    #if defined(dRenderVariant_marking)
        let thicknessAlpha = material.a < 1.0;
    #else
        let thicknessAlpha = uGlobal.uRenderMask == MaskTransparent;
    #endif
    if (thicknessAlpha && uObject.uAlphaThickness > 0.0) {
        material.a *= min(1.0, in.vRadius / (uObject.uAlphaThickness * uGlobal.uModelScale));
    }

    ${check_transparency_wgsl}

    #if defined(dRenderVariant_pick)
        ${pick_output_wgsl}
    #elif defined(dRenderVariant_marking)
        ${marking_output_wgsl}
    #elif defined(dRenderVariant_emissive)
        ${emissive_output_wgsl}
    #elif defined(dRenderVariant_depth)
        ${depth_output_wgsl()}
    #else
        ${apply_interior_wgsl}
        ${apply_light_wgsl}
        ${apply_marker_wgsl}
        #if defined(dRenderVariant_tracing)
            ${tracing_output_wgsl('normal')}
        #else
            ${apply_fog_wgsl}
            ${wboit_output_wgsl}
            ${dpoit_output_wgsl}
            out.color = color;
        #endif
    #endif
    return out;
}
`;
