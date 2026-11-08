/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { animation_wgsl, assign_wiggle_wgsl } from './chunks/animation.wgsl.js';
import { clip_flag_wgsl, clip_functions_wgsl, clip_instance_wgsl, clip_pixel_wgsl } from './chunks/clip.wgsl.js';
import { segment_wgsl } from './chunks/segment.wgsl.js';
import { sphere_shared_wgsl } from './chunks/sphere.wgsl.js';
import { common_wgsl } from './chunks/common.wgsl.js';
import { light_wgsl } from './chunks/light.wgsl.js';
import {
  apply_interior_wgsl,
  apply_light_wgsl,
  apply_fog_wgsl,
  apply_marker_wgsl,
  assign_material_wgsl,
  apply_occlusion_wgsl,
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
  assign_lod_factor_wgsl,
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
${sphere_shared_wgsl}

#if !defined(dClipPrimitive) && defined(dClipVariant_pixel) && defined(dClipObjects)
    #define dClipSphere
#endif

${vertex_output_wgsl(`    // the sphere's view-space center (negated) and radius, the same at every corner
    @location(9) @interpolate(flat) vSphere: vec4<f32>,
    @location(10) vPoint: vec3<f32>,
#if defined(dClipSphere)
    // 1 where the clip objects may cut the sphere (\`clipSphereClass\`): the fragment stage's clip tests skipped if 0
    @location(11) @interpolate(flat) vClipSome: u32,
#endif`)}

/**
 * The quad corner \`mapping\` of the bounding rectangle of a perspective-projected 3D sphere (\`sphereBounds\`), with the
 * off-axis terms of asymmetric projections (stereo, XR), \`uProjection[2].xy\`, zero for symmetric ones.
 */
fn sphereProjection(p: vec3<f32>, r: f32, mapping: vec2<f32>, w: f32) -> vec2<f32> {
    let projection = uGlobal.uProjection;
    let b = sphereBounds(p, r, projection);
    let xy = vec2<f32>(b.z + b.x, b.w + b.y) * -0.5 - mapping * vec2<f32>(b.z - b.x, b.w - b.y) * 0.5;
    return (xy - projection[2].xy) * w;
}

@vertex
fn vertexMain(in: VertexInput) -> VertexOutput {
    var out: VertexOutput;

    #if defined(dSphereResolve)
        // one triangle covering the viewport: its pixels shaded from the visibility buffer (\`resolveSphere\`)
        let corner = vec2<f32>(f32((in.vertexIndex << 1u) & 2u), f32(in.vertexIndex & 2u));
        out.position = vec4<f32>(corner * 2.0 - 1.0, 0.0, 1.0);
        return out;
    #else

    ${assign_quad_position_group_wgsl}
    ${assign_wiggle_wgsl}
    position = applyWiggle(position, group, wiggle);
    let aTransform = affineTransform(in.aTransform0, in.aTransform1, in.aTransform2);
    let transform = uObject.uObjectTransform * applyTumble(aTransform, invariantBoundingSphere(in));

    ${assign_varyings_wgsl}
    ${assign_size_wgsl}

    let modelScale = uGlobal.uModelScale;
    var radius = size * modelScale;

    let position4 = vec4<f32>(position, 1.0);
    out.vModelPosition = (uGlobal.uModel * transform * position4).xyz; // for clipping in the fragment stage

    ${assign_lod_factor_wgsl}
    radius *= lodFactor;

    var mvPosition = uGlobal.uModelView * transform * position4;
    out.vViewPosition = mvPosition.xyz;

    // through the near plane with none of it beyond
    var nearNone = false;
    #if defined(dApproximate)
        out.position = uGlobal.uProjection * vec4<f32>(mvPosition.xy + mapping * radius, mvPosition.z, 1.0);
    #else
        if (uGlobal.uIsOrtho == 1.0) {
            out.position = uGlobal.uProjection * vec4<f32>(mvPosition.xy + mapping * radius, mvPosition.z, 1.0);
        } else if (-mvPosition.z - radius > uGlobal.uNear) {
            let p = uGlobal.uProjection * vec4<f32>(mvPosition.xyz, 1.0);
            out.position = vec4<f32>(sphereProjection(mvPosition.xyz, radius, mapping, p.w), p.zw);
        } else {
            // through the near plane: the bounds of its part beyond it (the whole sphere's are unbounded where it
            // reaches behind the eye, a quad over the view), at a depth within the clip volume (the fragments write
            // theirs; the corners on one view plane, so the rays between them interpolate exactly)
            let b = sphereNearBounds(mvPosition.xyz, radius, uGlobal.uNear, uGlobal.uProjection);
            let xy = vec2<f32>(b.z + b.x, b.w + b.y) * -0.5 - mapping * vec2<f32>(b.z - b.x, b.w - b.y) * 0.5;
            out.position = vec4<f32>(xy - uGlobal.uProjection[2].xy, 0.5, 1.0);
            nearNone = b.x > b.z;
        }
    #endif

    let point4 = uGlobal.uInvProjection * out.position;
    out.vPoint = point4.xyz / point4.w;
    out.vSphere = vec4<f32>(-mvPosition.xyz / mvPosition.w, radius);

    // in front of the near plane (clip z is reversed)
    if (out.position.z > out.position.w) {
        mvPosition.z -= 2.0 * radius; // avoid clipping
        out.position.z = (uGlobal.uProjection * vec4<f32>(mvPosition.xyz, 1.0)).z;
    }

    if (lodCulled || nearNone) {
        // move out of [0, w] to 'discard' in the vertex stage
        out.position.z = 2.0 * out.position.w;
    }

    #if defined(dRenderVariant_marking)
        // the marking depth draws unmarked spheres only (\`uMarkingType\` 1), the mask marked ones: the others moved out of
        // [0, w], not shaded to be discarded
        let vertexMarker = select(floor(out.vMarkerClipping.x * 255.0 + 0.5), uObject.uMarker, uObject.uMarker != -1.0);
        if ((uGlobal.uMarkingType == 1) == (vertexMarker > 0.0)) {
            out.position.z = 2.0 * out.position.w;
        }
    #endif

    #if defined(dRenderVariant_colorBlended) || defined(dRenderVariant_colorWboit) || defined(dRenderVariant_colorDpoit)
        // drawn by the sphere raster (\`rasterizesSphere\`), with a margin: a sphere at its limits drawn by both, never
        // by neither (the raster's arithmetic is not the vertex stage's)
        if (uObject.uRasterMaxRadius > 0.0 && uGlobal.uIsOrtho == 0.0) {
            let b = spherePixelBounds(out.vViewPosition, radius, uGlobal.uProjection, uGlobal.uViewport.zw);
            if (rasterizesSphere(out.vViewPosition, radius, uGlobal.uNear * 1.001, b, uObject.uRasterMaxRadius * 0.999)) {
                out.position.z = 2.0 * out.position.w;
            }
        }
    #endif

    #if defined(dClipPrimitive) && !defined(dClipVariant_instance) && defined(dClipObjects)
        // the ghost pass (\`uClipGhost\`) keeps what the clip objects remove
        if (clipTest(out.vModelPosition / modelScale, clipFlag) != (uGlobal.uClipGhost != 0u)) {
            // move out of [0, w] to 'discard' in the vertex stage
            out.position.z = 2.0 * out.position.w;
        }
    #elif defined(dClipSphere)
        // pixel clipping: a sphere clipped whole moved out of [0, w], one not cut not tested per fragment; the ghost
        // pass (\`uClipGhost\`): a sphere whole by its center, kept where clipped, not tested per fragment
        var clipClass: u32;
        if (uGlobal.uClipGhost != 0u) {
            clipClass = select(2u, 0u, clipTest(out.vModelPosition / modelScale, clipFlag));
        } else {
            clipClass = clipSphereClass(out.vModelPosition / modelScale, radius / modelScale, clipFlag);
        }
        if (clipClass == 2u) {
            out.position.z = 2.0 * out.position.w;
        }
        out.vClipSome = select(0u, 1u, clipClass == 1u);
    #else
        ${clip_instance_wgsl}
    #endif
    return out;
    #endif
}

#if defined(dSphereResolve)
struct Resolved {
    out: VertexOutput,
    hit: bool,
};

/**
 * The sphere of the pixel at \`position\` in the visibility buffer as the vertex stage outputs it for the impostor
 * (\`hit\` false where none of the object's): its chunk's instance and item (\`SphereRasterTable\`), sized and levelled
 * as the raster drew it, the theme layers read by its indices (\`assign_varyings\`)
 */
fn resolveSphere(position: vec4<f32>) -> Resolved {
    var resolved: Resolved;
    resolved.hit = false;
    let pixel = vec2<u32>(position.xy);
    let v = sRasterVisibility[pixel.y * u32(uObject.uRasterWidth) + pixel.x];
    if (v.x == 0u) {
        return resolved;
    }
    let chunk = (v.x - 1u) / 64u;
    let base = u32(uObject.uRasterChunkBase);
    if (chunk < base || chunk >= base + u32(uObject.uRasterChunkCount)) {
        return resolved;
    }
    let e = sRasterChunks[chunk * 2u];
    let instance = sRasterChunks[chunk * 2u + 1u];
    let vertexId = i32(e.y + (v.x - 1u) % 64u);
    let positionGroup = read_tPositionGroup(vertexId);
    let group = positionGroup.w;
    let o = e.x * 3u;
    // the object's transient transform after the instance's, as the raster applies it
    let transform = uObject.uObjectTransform * affineTransform(sRasterTransforms[o], sRasterTransforms[o + 1u], sRasterTransforms[o + 2u]);

    // the instance as the vertex stage reads it: the group offsets the chunk holds
    var in: VertexInput;
    in.aInstance = f32(instance.x);
    #if defined(dSegmented)
        in.aSegment = vec4<f32>(0.0, f32(instance.y), f32(instance.z), 0.0);
    #endif
    var out: VertexOutput;
    ${assign_varyings_wgsl}
    ${assign_size_wgsl}

    let position4 = vec4<f32>(positionGroup.xyz, 1.0);
    out.vModelPosition = (uGlobal.uModel * transform * position4).xyz;
    var radius = size * uGlobal.uModelScale;
    // the chunk's level at its scale, as \`assign_lod_factor\`
    let lod = uObject.uRasterLevels[e.w];
    let scale = bitcast<f32>(instance.w);
    if (lod.w != 0.0 && (lod.x != 0.0 || lod.y != 0.0)) {
        if (uGlobal.uModelScale != 1.0) {
            radius *= scale;
        } else {
            let d = (dot(uGlobal.uCameraPlane.xyz, out.vModelPosition) + uGlobal.uCameraPlane.w) / uGlobal.uModelScale;
            radius *= lodAt(d, lod.xyz, scale).factor;
        }
    }
    // the center as the raster transforms it: the instance's transform first
    let mvPosition = uGlobal.uModelView * (transform * position4);
    out.vViewPosition = mvPosition.xyz;
    out.vSphere = vec4<f32>(-mvPosition.xyz / mvPosition.w, radius);
    // the view ray through the pixel's center as the raster casts it (\`sphere-raster.wgsl.ts\`), as the quad's point
    let viewport = uGlobal.uViewport.zw;
    let ndc = vec2<f32>(position.x / viewport.x * 2.0 - 1.0, 1.0 - position.y / viewport.y * 2.0);
    let projection = uGlobal.uProjection;
    out.vPoint = vec3<f32>(
        (ndc.x + projection[2].x) / projection[0][0],
        (ndc.y + projection[2].y) / projection[1][1],
        -1.0
    );
    out.position = position;
    #if defined(dClipSphere)
        // the raster tested the sphere: every one by pixel, as it drew it
        out.vClipSome = 1u;
    #endif
    resolved.out = out;
    resolved.hit = true;
    return resolved;
}
#endif

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

/** \`clipSome\`: the clip objects may cut the sphere (\`vClipSome\`), else its clip tests are skipped */
fn sphereImpostor(radius: f32, point: vec3<f32>, pointViewPosition: vec3<f32>, clipFlag: i32, clipSome: bool) -> SphereHit {
    var h: SphereHit;
    h.hit = false;
    h.isCap = false;

    let isOrtho = uGlobal.uIsOrtho;
    let cameraSpherePos = -pointViewPosition;
    let rayOrigin = mix(vec3<f32>(0.0), point, isOrtho);
    let rayDirection = mix(normalize(point), vec3<f32>(0.0, 0.0, 1.0), isOrtho);
    let cameraSphereDir = mix(cameraSpherePos, rayOrigin - cameraSpherePos, isOrtho);

    let hits = raySphere(rayDirection, cameraSphereDir, radius);
    // the resolve's pixels are the raster's hits: a grazing miss here is its rounding, the tangent point is taken
    #if !defined(dSphereResolve)
        if (hits.x < 0.0) {
            return h;
        }
    #endif
    let posT = mix(hits.z, hits.y, isOrtho);
    let negT = mix(hits.y, hits.z, isOrtho);

    #if defined(dRenderVariant_depth)
        if (uGlobal.uDepthBack != 0u) {
            // back-depth pass wants the farthest point that is not clipped
            h.viewPosition = rayDirection * posT + rayOrigin;
            h.modelPosition = (uGlobal.uInvView * vec4<f32>(h.viewPosition, 1.0)).xyz;
            h.normal = -normalize(h.viewPosition - cameraSpherePos);
            h.interior = true;
            #if defined(dSolidInterior) && defined(dClipSphere)
            if (clipSome) {
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
        if (clipSome && clipTest(h.modelPosition / uGlobal.uModelScale, clipFlag)) {
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
                    var sCap = sNear;
                    if (clipSome) {
                        sCap = clipExit(frontModelPosition / uGlobal.uModelScale, h.modelPosition / uGlobal.uModelScale, sNear, clipFlag);
                    }
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
fn fragmentMain(fragmentIn: VertexOutput) -> FragmentOutput {
    #if defined(dSphereResolve)
        let resolved = resolveSphere(fragmentIn.position);
        if (!resolved.hit) {
            discard;
        }
        let in = resolved.out;
    #else
        let in = fragmentIn;
    #endif
    var out: FragmentOutput;
    ${clip_flag_wgsl('in.vMarkerClipping.y')}

    #if defined(dApproximate)
        var pointDir = -in.vSphere.xyz - in.vPoint;
        if (dot(pointDir, pointDir) > in.vSphere.w * in.vSphere.w) {
            discard;
        }
        let viewPosition = -in.vSphere.xyz;
        let modelPosition = in.vModelPosition;
        let fragmentDepth = in.position.z;
        pointDir.z -= cos(length(pointDir)) * in.vSphere.w * 0.5;
        let cameraNormal = -normalize(pointDir);
        let interior = false;
        let isCap = false;
    #else
        #if defined(dClipSphere)
            let clipSome = in.vClipSome != 0u;
        #else
            let clipSome = false;
        #endif
        let h = sphereImpostor(in.vSphere.w, in.vPoint, in.vSphere.xyz, clipFlag, clipSome);
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
    #if defined(dClipSphere)
        if (!isCap && in.vClipSome != 0u) {
            ${clip_pixel_wgsl('modelPosition')}
        }
    #elif !defined(dClipPrimitive)
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
        material.a *= min(1.0, in.vSphere.w / (uObject.uAlphaThickness * uGlobal.uModelScale));
    }

    ${apply_occlusion_wgsl}
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
