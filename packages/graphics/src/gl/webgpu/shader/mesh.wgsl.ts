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
import { iso_vertex_wgsl } from './marching-cubes-chunks.wgsl.js';
import { packing_wgsl } from './chunks/packing.wgsl.js';
import { read_from_texture_wgsl } from './chunks/read-from-texture.wgsl.js';
import { assign_varyings_wgsl, grid_functions_wgsl, vertex_output_wgsl } from './chunks/varyings.wgsl.js';

export const mesh_wgsl = /* wgsl */ `
${common_wgsl}
${read_from_texture_wgsl}
${animation_wgsl}
${light_wgsl}
${material_functions_wgsl}
${segment_wgsl}
${clip_functions_wgsl}
${packing_wgsl}
${grid_functions_wgsl}

${vertex_output_wgsl('    @location(9) vNormal: vec3<f32>,')}
${fragment_output_wgsl(`#if defined(dSolidInterior)
    @builtin(frag_depth) fragDepth: f32,
#endif`)}

#if defined(dGeometryType_textureMesh) && defined(dIsoChunks)
${iso_vertex_wgsl('sIsoTable', 'sIsoVertices')}
#endif

@vertex
fn vertexMain(in: VertexInput) -> VertexOutput {
    var out: VertexOutput;

    let vertexId = i32(in.vertexIndex);
    #if defined(dGeometryType_textureMesh) && defined(dIsoChunks)
        // unindexed, the vertex data in pages of the vertex arena; vertices beyond the drawn count (the page table's
        // fourth word, e.g. of draws of the CPU's upper bound) are left out below
        let isDrawn = in.vertexIndex < isoDrawn();
        let slot = select(0u, isoSlot(in.vertexIndex), isDrawn);
        let group = f32(isoPackedGroup(slot)) - 1.0;
        let basePosition = isoPosition(slot);
        let baseNormal = isoNormal(slot);
    #elif defined(dGeometryType_textureMesh)
        // unindexed, the vertex data in textures (e.g. computed on the GPU)
        let group = unpackRGBToInt(readFromTexture(tGroup, vertexId).rgb);
        let basePosition = readFromTexture(tPosition, vertexId).xyz;
        let baseNormal = readFromTexture(tNormal, vertexId).xyz;
    #else
        let group = in.aGroup;
        let basePosition = in.aPosition;
        let baseNormal = in.aNormal;
    #endif

    let aTransform = mat4x4<f32>(in.aTransform0, in.aTransform1, in.aTransform2, in.aTransform3);
    let transform = applyTumble(aTransform, in.aInstance, f32(uObject.uObjectId), invariantBoundingSphere(in));
    let model = uGlobal.uModel * transform;
    let modelView = uGlobal.uView * model;
    ${assign_wiggle_wgsl}
    let position4 = vec4<f32>(applyWiggle(basePosition, group, wiggle), 1.0);
    out.vModelPosition = (model * position4).xyz;
    let mvPosition = modelView * position4;
    out.vViewPosition = mvPosition.xyz;
    out.position = uGlobal.uProjection * mvPosition;
    #if defined(dGeometryType_textureMesh) && defined(dIsoChunks)
        if (!isDrawn) {
            // outside the clip volume: whole triangles (the drawn count is a multiple of 3) are clipped
            out.position = vec4<f32>(0.0, 0.0, 2.0, 1.0);
        }
    #endif

    // after the positions: grid data (smoothed colors) is read at the vertex position
    let gridPosition = basePosition;
    ${assign_varyings_wgsl}

    var normal = normalize(adjoint(modelView) * normalize(baseNormal));
    #if defined(dFlipSided)
        if (uMaterial.uDoubleSided == 0u) {
            normal = -normal;
        }
    #endif
    out.vNormal = normal;

    #if defined(dSolidInterior)
        if (uDraw.uSolidInteriorPass != 0) {
            // keep the back faces of the cap passes in front of the far plane
            // in front of the far plane (clip z is reversed)
            out.position.z = max(out.position.z, out.position.w * 0.0001);
        }
    #endif
    ${clip_instance_wgsl}
    return out;
}

@fragment
fn fragmentMain(in: VertexOutput, @builtin(front_facing) frontFacing: bool) -> FragmentOutput {
    var out: FragmentOutput;
    ${clip_flag_wgsl('in.vClipping')}

    // derivatives in uniform control flow; framebuffer y points down, so dpdy has the opposite sign of GL's dFdy
    let fdx = dpdx(in.vViewPosition);
    let fdy = dpdy(in.vViewPosition);

    #if defined(dSolidInterior)
        // solid interior caps: the back faces seen through the cut (near plane or clip object), drawn where the stencil
        // of the mark pass (\`uSolidInteriorPass\` 2) counts more back than front faces, at the cut (pass 1)
        var fragmentDepth = in.position.z;
        let capPass = uDraw.uSolidInteriorPass != 0;
        var viewPosition = in.vViewPosition;
        var modelPosition = in.vModelPosition;
        if (capPass) {
            let nearZ = -uGlobal.uNear * 1.0001;
            let vp = in.vViewPosition;
            let nearPosition = vec3<f32>(mix(vp.xy * (nearZ / vp.z), vp.xy, uGlobal.uIsOrtho), nearZ);
            var s = 0.0;
            #if defined(dClipObjects)
                if (uDraw.uSolidInteriorClip >= 0) {
                    // exempt parts are not cut by this clip object, so they neither count nor fill
                    if (clipIgnored(uDraw.uSolidInteriorClip, clipFlag)) {
                        discard;
                    }
                    let nearModelPosition = (uGlobal.uInvView * vec4<f32>(nearPosition, 1.0)).xyz;
                    s = clipCapExit(nearModelPosition / uGlobal.uModelScale, in.vModelPosition / uGlobal.uModelScale, uGlobal.uDepthBack != 0u);
                    if (s < 0.0 || s > 1.0) {
                        discard;
                    }
                }
            #endif
            if (uDraw.uSolidInteriorPass == 2) {
                // stencil marks only
                out.fragDepth = fragmentDepth;
                return out;
            }
            viewPosition = mix(nearPosition, in.vViewPosition, s);
            modelPosition = (uGlobal.uInvView * vec4<f32>(viewPosition, 1.0)).xyz;
            fragmentDepth = mix(calcDepth(viewPosition), in.position.z, 0.0001);
            // beyond the far plane (depth is reversed)
            if (fragmentDepth < 0.0) {
                discard;
            }
        }
        out.fragDepth = fragmentDepth;
        ${clip_pixel_wgsl('modelPosition')}
        // after the mark pass returned (it needs every face) and at the cap position
        if (fadeLod(in.position, modelPosition)) {
            discard;
        }
    #else
        if (fadeLod(in.position, in.vModelPosition)) {
            discard;
        }
        ${clip_pixel_wgsl('in.vModelPosition')}

        let fragmentDepth = in.position.z;
        let modelPosition = in.vModelPosition;
        let viewPosition = in.vViewPosition;
    #endif

    var interior = !frontFacing;
    #if defined(dSolidInterior)
        // the fill may draw front faces (reflected instances), a cap is always interior
        #if defined(dRenderVariant_depth)
            // a cap is the visible surface, so the depth passes must not discard it as interior
            if (capPass) {
                interior = false;
            }
        #else
            if (capPass) {
                interior = true;
            }
        #endif
    #endif

    #if defined(dFlatShaded)
        var normal = normalize(cross(fdx, fdy));
    #else
        var normal = -normalize(in.vNormal);
        if (uMaterial.uDoubleSided != 0u) {
            normal *= select(-1.0, 1.0, frontFacing);
        }
    #endif
    #if defined(dFlipSided)
        normal = -normal;
    #endif
    #if defined(dSolidInterior)
        if (capPass) {
            normal = vec3<f32>(0.0, 0.0, -1.0);
            #if defined(dClipObjects)
                if (uDraw.uSolidInteriorClip >= 0) {
                    let invView = uGlobal.uInvView;
                    let invView3 = mat3x3<f32>(invView[0].xyz, invView[1].xyz, invView[2].xyz);
                    normal = normalize(clipCapNormal(modelPosition / uGlobal.uModelScale) * invView3);
                }
            #endif
        }
    #endif

    ${assign_material_wgsl}
    ${check_transparency_wgsl}

    #if defined(dRenderVariant_pick)
        ${pick_output_wgsl}
        #if defined(dSolidInterior)
            if (capPass) {
                // 16777214 is PickingId.Null, the cap picks the whole instance
                out.group = vec4<f32>(packIntToRGB(16777214.0), 1.0);
            }
        #endif
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
