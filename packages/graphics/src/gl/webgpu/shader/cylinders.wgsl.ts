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
import { assign_size_wgsl, assign_varyings_wgsl, vertex_output_wgsl } from './chunks/varyings.wgsl.js';

/** Cylinder impostors: a box per cylinder (6 vertices, 4 triangles), ray-cylinder intersection per fragment */
export const cylinders_wgsl = /* wgsl */ `
${common_wgsl}
${read_from_texture_wgsl}
${animation_wgsl}
${light_wgsl}
${material_functions_wgsl}
${segment_wgsl}
${clip_functions_wgsl}
${packing_wgsl}

${vertex_output_wgsl(`    @location(9) vStart: vec3<f32>,
    @location(10) vEnd: vec3<f32>,
    @location(11) vSize: f32,
    @location(12) vCap: f32,`)}

${fragment_output_wgsl('    @builtin(frag_depth) fragDepth: f32,')}

@vertex
fn vertexMain(in: VertexInput) -> VertexOutput {
    var out: VertexOutput;

    let vertexId = i32(in.vertexIndex);
    let group = in.aGroup;
    ${assign_varyings_wgsl}
    ${assign_size_wgsl}

    let aTransform = mat4x4<f32>(in.aTransform0, in.aTransform1, in.aTransform2, in.aTransform3);
    let transform = applyTumble(aTransform, in.aInstance, f32(uObject.uObjectId), invariantBoundingSphere(in));
    ${assign_wiggle_wgsl}
    let wigStart = applyWiggle(in.aStart, group, wiggle);
    let wigEnd = applyWiggle(in.aEnd, group, wiggle);
    let modelTransform = uGlobal.uModel * transform;

    let start = (modelTransform * vec4<f32>(wigStart, 1.0)).xyz;
    let end = (modelTransform * vec4<f32>(wigEnd, 1.0)).xyz;
    let vSize = size * in.aScale * uGlobal.uModelScale;
    out.vStart = start;
    out.vEnd = end;
    out.vSize = vSize;
    out.vCap = in.aCap;

    var modelPosition = (start + end) * 0.5;
    let camDir = -mix(normalize(modelPosition - uGlobal.uCameraPosition), uGlobal.uCameraDir, uGlobal.uIsOrtho);
    var dir = end - start;
    var f = select(0.0, 1.0, in.aMapping.x > 0.0);
    // ensure cylinder 'dir' is pointing towards the camera
    if (dot(camDir, dir) < 0.0) {
        dir = -dir;
        f = 1.0 - f;
    }

    var left = cross(camDir, dir);
    var up = cross(left, dir);
    left = vSize * normalize(left);
    up = vSize * normalize(up);

    // move vertex in object-space from center to corner
    modelPosition += in.aMapping.x * dir + in.aMapping.y * left + in.aMapping.z * up;
    out.vModelPosition = modelPosition;

    let mvPosition = uGlobal.uView * vec4<f32>(modelPosition, 1.0);
    out.vViewPosition = mvPosition.xyz;
    out.position = uGlobal.uProjection * mvPosition;

    // in front of the near plane (clip z is reversed): onto it, to avoid clipping the box; the fragment depth is the
    // impostor's
    if (out.position.z > out.position.w) {
        out.position.z = out.position.w;
    }

    #if defined(dDualColor) && (defined(dRenderVariant_color) || defined(dRenderVariant_tracing) || defined(dRenderVariant_emissive)) && defined(dColorType_texture)
        // dual-color mixing (group color types)
        // - for aColorMode between 0 and 1 use aColorMode to interpolate
        // - for aColorMode == 2 do nothing, i.e., use vColor
        // - for aColorMode == 3 use position on cylinder axis to interpolate
        if (dualColor) {
            if (in.aColorMode <= 1.0) {
                out.vColor = vec4<f32>(mix(out.vColor.rgb, color2, in.aColorMode), out.vColor.a);
            } else if (in.aColorMode == 3.0) {
                out.vColor = vec4<f32>(mix(out.vColor.rgb, color2, mix(-0.25, 1.25, f / 1.5)), out.vColor.a);
            }
        }
    #endif
    ${clip_instance_wgsl}
    return out;
}

struct CylinderHit {
    hit: bool,
    cameraNormal: vec3<f32>,
    interior: bool,
    modelPosition: vec3<f32>,
    viewPosition: vec3<f32>,
    depth: f32,
    /** on a clip cap of a solid interior, not tested by pixel clipping */
    isCap: bool,
};

#if defined(dClipVariant_pixel) && defined(dClipObjects)
    #define dClipCylinder
#endif

fn cylinderHitAt(rayOrigin: vec3<f32>, rayDir: vec3<f32>, t: f32, cameraNormal: vec3<f32>, interior: bool) -> CylinderHit {
    var h: CylinderHit;
    h.cameraNormal = cameraNormal;
    h.interior = interior;
    h.modelPosition = rayOrigin + t * rayDir;
    h.viewPosition = (uGlobal.uView * vec4<f32>(h.modelPosition, 1.0)).xyz;
    h.depth = calcDepth(h.viewPosition);
    h.hit = h.depth < 1.0;
    h.isCap = false;
    return h;
}

/** Where the ray enters the cylinder, for the caps of solid interiors */
struct CylinderFront {
    t: f32,
    depth: f32,
    modelPosition: vec3<f32>,
    /** removed by a clip object */
    clipped: bool,
};

/** Whether a clip object removes \`modelPosition\` of the front */
fn isFrontClipped(modelPosition: vec3<f32>, clipFlag: i32) -> bool {
    #if defined(dClipCylinder)
        return clipTest(modelPosition / uGlobal.uModelScale, clipFlag);
    #else
        return false;
    #endif
}

#if defined(dSolidInterior)
    /**
     * The interior hit, or where the front is cut away (by the near plane or a clip object), the cap there; a miss when
     * clip objects remove the whole segment
     */
    fn solidInteriorCap(hit: CylinderHit, rayOrigin: vec3<f32>, rayDir: vec3<f32>, t: f32, front: CylinderFront, radius: f32, clipFlag: i32) -> CylinderHit {
        var h = hit;
        let cameraRayOrigin = (uGlobal.uView * vec4<f32>(rayOrigin, 1.0)).xyz;
        let cameraRayDir = (uGlobal.uView * vec4<f32>(rayDir, 0.0)).xyz;
        let nearT = -(uGlobal.uNear + cameraRayOrigin.z) / cameraRayDir.z;
        let sNear = select(clamp((nearT - front.t) / (t - front.t), 0.0, 1.0), 0.0, front.depth < 1.0);
        #if defined(dClipCylinder)
            var sCap = sNear;
            if (front.clipped || front.depth >= 1.0) {
                sCap = clipExit(front.modelPosition / uGlobal.uModelScale, h.modelPosition / uGlobal.uModelScale, sNear, clipFlag);
            }
            if (sCap < 0.0) {
                h.hit = false;
                return h;
            }
        #else
            let sCap = sNear;
        #endif
        if (sCap == sNear && front.depth >= 1.0) {
            h.depth = 1.0 - (0.0000002 / radius);
            h.cameraNormal = -rayDir;
            h.viewPosition = cameraRayOrigin + nearT * cameraRayDir;
            h.modelPosition = (uGlobal.uInvView * vec4<f32>(h.viewPosition, 1.0)).xyz;
        }
        #if defined(dClipCylinder)
            else if (front.clipped || front.depth >= 1.0) {
                h.modelPosition = mix(front.modelPosition, h.modelPosition, sCap);
                h.viewPosition = (uGlobal.uView * vec4<f32>(h.modelPosition, 1.0)).xyz;
                h.depth = calcDepth(h.viewPosition) - (0.0000002 / radius);
                h.cameraNormal = -clipNormal(h.modelPosition / uGlobal.uModelScale, clipFlag);
                h.isCap = true;
            }
        #endif
        return h;
    }

    /** The back-depth pass's far hit, its exit point moved towards the entry until a clip object no longer removes it */
    fn solidInteriorBackCap(hit: CylinderHit, entryPosition: vec3<f32>, clipFlag: i32) -> CylinderHit {
        var h = hit;
        #if defined(dClipCylinder)
            let s = clipExit(h.modelPosition / uGlobal.uModelScale, entryPosition / uGlobal.uModelScale, 0.0, clipFlag);
            if (s < 0.0) {
                h.hit = false;
                return h;
            }
            if (s > 0.0) {
                h.modelPosition = mix(h.modelPosition, entryPosition, s);
                h.viewPosition = (uGlobal.uView * vec4<f32>(h.modelPosition, 1.0)).xyz;
                h.depth = calcDepth(h.viewPosition);
                h.cameraNormal = -clipNormal(h.modelPosition / uGlobal.uModelScale, clipFlag);
                h.isCap = true;
            }
        #endif
        h.hit = h.depth < 1.0;
        return h;
    }
#else
    fn solidInteriorCap(hit: CylinderHit, rayOrigin: vec3<f32>, rayDir: vec3<f32>, t: f32, front: CylinderFront, radius: f32, clipFlag: i32) -> CylinderHit {
        return hit;
    }

    fn solidInteriorBackCap(hit: CylinderHit, entryPosition: vec3<f32>, clipFlag: i32) -> CylinderHit {
        return hit;
    }
#endif

// adapted from https://www.shadertoy.com/view/4lcSRn
// The MIT License, Copyright 2016 Inigo Quilez
fn cylinderImpostor(rayOrigin: vec3<f32>, rayDir: vec3<f32>, start: vec3<f32>, end: vec3<f32>, radius: f32, cap: f32, clipFlag: i32) -> CylinderHit {
    var r: CylinderHit;
    r.hit = false;

    let ba = end - start;
    let oc = rayOrigin - start;

    let baba = dot(ba, ba);
    let bard = dot(ba, rayDir);
    let baoc = dot(ba, oc);

    let k2 = baba - bard * bard;
    let k1 = baba * dot(oc, rayDir) - baoc * bard;
    let k0 = baba * dot(oc, oc) - baoc * baoc - radius * radius * baba;

    var h = k1 * k1 - k2 * k0;
    if (h < 0.0) {
        return r;
    }

    var topCap = (cap > 0.9 && cap < 1.1) || cap >= 2.9;
    var bottomCap = (cap > 1.9 && cap < 2.1) || cap >= 2.9;

    #if defined(dSolidInterior)
        let topInterior = !topCap;
        let bottomInterior = !bottomCap;
        topCap = true;
        bottomCap = true;
        let solidInterior = true;
    #else
        let topInterior = false;
        let bottomInterior = false;
        let solidInterior = false;
    #endif

    var clipped = false;
    // no front yet: at the far plane (depth is reversed)
    var front = CylinderFront(0.0, 0.0, rayOrigin, false);

    #if defined(dRenderVariant_depth)
        // the back-depth pass wants the farthest intersection
        let depthBack = uGlobal.uDepthBack != 0u;
    #else
        let depthBack = false;
    #endif

    // body outside
    h = sqrt(h);
    var t = (-k1 - h) / k2;
    var y = baoc + t * bard;
    // entry into the closed cylinder: the body entry clamped to the end-cap slab (for back caps of solid interiors)
    let entryT = select(t, max(t, min(-baoc / bard, (baba - baoc) / bard)), abs(bard) > 0.0);
    let entryPosition = rayOrigin + entryT * rayDir;
    if (!depthBack && y > 0.0 && y < baba) {
        r = cylinderHitAt(rayOrigin, rayDir, t, (oc + t * rayDir - ba * y / baba) / radius, false);
        front = CylinderFront(t, r.depth, r.modelPosition, isFrontClipped(r.modelPosition, clipFlag));
        if (front.clipped) {
            r.depth = 2.0;
            r.hit = false;
        }
        if (r.hit) {
            return r;
        }
        clipped = true;
    }

    if (!clipped && !depthBack) {
        var capT = 0.0;
        var capNormal = vec3<f32>(0.0);
        var capInterior = false;
        var hasCap = false;
        if (topCap && y < 0.0) {
            // top cap
            capT = -baoc / bard;
            capNormal = -ba / baba;
            capInterior = topInterior;
            hasCap = abs(k1 + k2 * capT) < h;
        } else if (bottomCap && y >= 0.0) {
            // bottom cap
            capT = (baba - baoc) / bard;
            capNormal = ba / baba;
            capInterior = bottomInterior;
            hasCap = abs(k1 + k2 * capT) < h;
        }
        if (hasCap) {
            r = cylinderHitAt(rayOrigin, rayDir, capT, capNormal, capInterior);
            front = CylinderFront(capT, r.depth, r.modelPosition, isFrontClipped(r.modelPosition, clipFlag));
            if (front.clipped) {
                r.depth = 2.0;
                r.hit = false;
            }
            if (r.hit) {
                #if defined(dSolidInterior)
                    if (r.interior) {
                        r.cameraNormal = -rayDir;
                    }
                #endif
                #if !defined(dClipCylinder)
                    // open ends of solid interiors are seen through (without clipping)
                    r.hit = !r.interior;
                #endif
                return r;
            }
        }
    }

    if (uMaterial.uDoubleSided != 0u || solidInterior || depthBack) {
        // body inside
        h = -h;
        t = (-k1 - h) / k2;
        y = baoc + t * bard;
        if (y > 0.0 && y < baba) {
            r = cylinderHitAt(rayOrigin, rayDir, t, -(oc + t * rayDir - ba * y / baba) / radius, true);
            if (r.hit) {
                if (depthBack) {
                    return solidInteriorBackCap(r, entryPosition, clipFlag);
                }
                return solidInteriorCap(r, rayOrigin, rayDir, t, front, radius, clipFlag);
            }
        }

        if (topCap && y < 0.0) {
            // top cap
            t = -baoc / bard;
            if (abs(k1 + k2 * t) < -h) {
                r = cylinderHitAt(rayOrigin, rayDir, t, ba / baba, true);
                if (r.hit) {
                    if (depthBack) {
                        return solidInteriorBackCap(r, entryPosition, clipFlag);
                    }
                    return solidInteriorCap(r, rayOrigin, rayDir, t, front, radius, clipFlag);
                }
            }
        } else if (bottomCap && y >= 0.0) {
            // bottom cap
            t = (baba - baoc) / bard;
            if (abs(k1 + k2 * t) < -h) {
                r = cylinderHitAt(rayOrigin, rayDir, t, -ba / baba, true);
                if (r.hit) {
                    if (depthBack) {
                        return solidInteriorBackCap(r, entryPosition, clipFlag);
                    }
                    return solidInteriorCap(r, rayOrigin, rayDir, t, front, radius, clipFlag);
                }
            }
        }
    }

    r.hit = false;
    return r;
}

@fragment
fn fragmentMain(in: VertexOutput) -> FragmentOutput {
    var out: FragmentOutput;
    ${clip_flag_wgsl('in.vClipping')}

    let rayOrigin = in.vModelPosition;
    let rayDir = mix(normalize(in.vModelPosition - uGlobal.uCameraPosition), uGlobal.uCameraDir, uGlobal.uIsOrtho);

    let h = cylinderImpostor(rayOrigin, rayDir, in.vStart, in.vEnd, in.vSize, in.vCap, clipFlag);
    if (!h.hit || h.depth < 0.0 || h.depth > 1.0) {
        discard;
    }
    out.fragDepth = h.depth;

    let modelPosition = h.modelPosition;
    let viewPosition = h.viewPosition;
    let fragmentDepth = h.depth;
    let interior = h.interior;

    if (fadeLod(in.position, modelPosition)) {
        discard;
    }
    if (!h.isCap) {
        ${clip_pixel_wgsl('modelPosition')}
    }

    var normal = normalize(adjoint(uGlobal.uView) * -normalize(h.cameraNormal));

    ${assign_material_wgsl}
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
