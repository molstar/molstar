/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { clip_functions_wgsl, clip_instance_wgsl } from './chunks/clip.wgsl.js';
import { segment_wgsl } from './chunks/segment.wgsl.js';
import { common_wgsl } from './chunks/common.wgsl.js';
import {
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
import { assign_varyings_wgsl, vertex_output_wgsl } from './chunks/varyings.wgsl.js';

/**
 * Gaussian splats as quads along the principal axes of each splat's projected 2D Gaussian (EWA): the 3D covariance
 * `R S Sᵀ Rᵀ` from the splat's rotation and scales (`tSplat`), projected by the Jacobian of the projection at its
 * center; unlit
 */
export const gaussian_splats_wgsl = /* wgsl */ `
${common_wgsl}
#if defined(dColorType_direct)
    // SOG's base color, decoded in the vertex stage (\`vColor\`)
    #define dColorType_varying
#endif
${read_from_texture_wgsl}
${material_functions_wgsl}
${segment_wgsl}
${clip_functions_wgsl}
${packing_wgsl}

${vertex_output_wgsl(`    @location(9) vQuad: vec2<f32>,
    @location(10) @interpolate(flat) vSplatOpacity: f32,
    @location(11) @interpolate(flat) vSplatSh: vec3<f32>,`)}
${fragment_output_wgsl()}

#if defined(dSplatEncoding_sog)
    /** a texel of a SOG image as bytes */
    fn sogBytes(v: vec4<f32>) -> vec4<f32> {
        return floor(v * 255.0 + 0.5);
    }

    /** splat \`i\`'s texel of SOG's images (a pixel a splat, row by row) */
    fn sogTexel(i: i32) -> vec2<i32> {
        let w = i32(uObject.uSogImageDim.x);
        return vec2<i32>(i % w, i / w);
    }

    /** the codebooks' entry \`b\`: scale (r), base color (g), view-dependent color (b) */
    fn sogCode(b: f32) -> vec4<f32> {
        return textureLoad(tSogCodebook, vec2<i32>(i32(b), 0), 0);
    }

    /** splat \`i\`'s center: 16 bits an axis within the means' range, out of the log domain */
    fn sogPosition(i: i32) -> vec3<f32> {
        let p = sogTexel(i);
        let lo = sogBytes(textureLoad(tSogMeansLow, p, 0)).rgb;
        let hi = sogBytes(textureLoad(tSogMeansHigh, p, 0)).rgb;
        let n = uObject.uSogMeansMin + (uObject.uSogMeansMax - uObject.uSogMeansMin) * ((lo + hi * 256.0) / 65535.0);
        return sign(n) * (exp(abs(n)) - 1.0);
    }

    /** splat \`i\`'s rotation (\`x, y, z, w\`): version 3 an exponential map, version 2 the smallest three */
    fn sogRotation(i: i32) -> vec4<f32> {
        let b = sogBytes(textureLoad(tSogQuats, sogTexel(i), 0));
        #if defined(dSogVersion3)
            let p = (b.xyz - 128.0) / 127.0;
            let r = length(p);
            let t = 1.5707963267948966 * r;
            let k = select(1.5707963267948966, sin(t) / r, r > 1e-12);
            return vec4<f32>(p * k, cos(t));
        #else
            // the largest component's index in \`w, x, y, z\` (252 to 255), the stored three the others in order
            let largest = b.w - 252.0;
            let c = (b.xyz / 255.0 * 2.0 - 1.0) * 0.7071067811865476;
            let m = sqrt(max(0.0, 1.0 - dot(c, c)));
            if (largest < -0.5 || largest > 3.5) {
                return vec4<f32>(0.0, 0.0, 0.0, 1.0);
            }
            if (largest < 0.5) {
                return vec4<f32>(c, m);
            }
            if (largest < 1.5) {
                return vec4<f32>(m, c.y, c.z, c.x);
            }
            if (largest < 2.5) {
                return vec4<f32>(c.y, m, c.z, c.x);
            }
            return vec4<f32>(c.y, c.z, m, c.x);
        #endif
    }

    /**
     * the palette's coefficient \`k\` of entry \`label\` (rgb, \`uSogShCoefficients\` pixels an entry, 64 entries a row),
     * see \`getSogShLayout\`; \`count\` unused (the palette's own)
     */
    fn splatSh(label: i32, count: i32, k: i32) -> vec3<f32> {
        #if defined(dSogSh)
            let u = (label % 64) * uObject.uSogShCoefficients + k;
            let b = sogBytes(textureLoad(tSogShCentroids, vec2<i32>(u, label / 64), 0)).rgb;
            return vec3<f32>(sogCode(b.r).b, sogCode(b.g).b, sogCode(b.b).b);
        #else
            return vec3<f32>(0.0);
        #endif
    }
#else
    /** splat \`splat\`'s coefficient of basis function \`k\` past the zeroth (rgb), see \`SplatSh\` */
    fn splatSh(splat: i32, count: i32, k: i32) -> vec3<f32> {
        return (read_tSplatSh(splat * count + k).rgb * 2.0 - 1.0) * uObject.uSplatShScale;
    }
#endif

/**
 * The view-dependent color of splat \`splat\` past its base color (3D Gaussian splatting's \`eval_sh\` without the zeroth),
 * \`dir\` the unit direction from the camera to it in its own frame
 */
fn evalSplatSh(splat: i32, dir: vec3<f32>) -> vec3<f32> {
    let degree = uObject.uSplatShDegree;
    let n = select(select(15, 8, degree == 2), 3, degree == 1);
    let x = dir.x;
    let y = dir.y;
    let z = dir.z;
    var c = 0.4886025119029199 * (-y * splatSh(splat, n, 0) + z * splatSh(splat, n, 1) - x * splatSh(splat, n, 2));
    if (degree > 1) {
        let xx = x * x;
        let yy = y * y;
        let zz = z * z;
        c += 1.0925484305920792 * x * y * splatSh(splat, n, 3)
            - 1.0925484305920792 * y * z * splatSh(splat, n, 4)
            + 0.31539156525252005 * (2.0 * zz - xx - yy) * splatSh(splat, n, 5)
            - 1.0925484305920792 * x * z * splatSh(splat, n, 6)
            + 0.5462742152960396 * (xx - yy) * splatSh(splat, n, 7);
        if (degree > 2) {
            c += -0.5900435899266435 * y * (3.0 * xx - yy) * splatSh(splat, n, 8)
                + 2.890611442640554 * x * y * z * splatSh(splat, n, 9)
                - 0.4570457994644658 * y * (4.0 * zz - xx - yy) * splatSh(splat, n, 10)
                + 0.3731763325901154 * z * (2.0 * zz - 3.0 * xx - 3.0 * yy) * splatSh(splat, n, 11)
                - 0.4570457994644658 * x * (4.0 * zz - xx - yy) * splatSh(splat, n, 12)
                + 1.445305721320277 * z * (xx - yy) * splatSh(splat, n, 13)
                - 0.5900435899266435 * x * (xx - 3.0 * yy) * splatSh(splat, n, 14);
        }
    }
    return c;
}

@vertex
fn vertexMain(in: VertexInput) -> VertexOutput {
    var out: VertexOutput;

    let m = i32(in.vertexIndex) % 6;
    var mapping = vec2<f32>(1.0, 1.0); // vertices 2 and 5
    if (m == 0) {
        mapping = vec2<f32>(-1.0, 1.0);
    } else if (m == 1 || m == 3) {
        mapping = vec2<f32>(-1.0, -1.0);
    } else if (m == 4) {
        mapping = vec2<f32>(1.0, -1.0);
    }
    // the splat of this draw slot: the instance's block of the order, back to front for the view (\`tSplatOrder\`, a
    // u32 in four bytes)
    let slot = i32(in.vertexIndex) / 6;
    #if defined(dSplatGpuOrder)
        // sorted on the GPU (\`sSplatOrder\`)
        let vertexId = i32(sSplatOrder[u32(i32(in.aInstance) * uObject.uVertexCount + slot)]);
    #else
        let orderBytes = vec4<i32>(read_tSplatOrder(i32(in.aInstance) * uObject.uVertexCount + slot) * 255.0 + 0.5);
        let vertexId = orderBytes.x + orderBytes.y * 256 + orderBytes.z * 65536 + orderBytes.w * 16777216;
    #endif
    #if defined(dSplatEncoding_sog)
        // SOG's images read in place (\`GaussianSplats.Sog\`): the group the splat's index
        let position = sogPosition(vertexId);
        let group = vertexId;
        let sh0 = sogBytes(textureLoad(tSogSh0, sogTexel(vertexId), 0));
        #if defined(dColorType_direct) && (defined(dRenderVariant_color) || defined(dRenderVariant_tracing) || defined(dRenderVariant_emissive))
            // the base color the file's (\`0.5 + SH_C0 dc\`), before overpaint mixes over it
            out.vColor = vec4<f32>(0.5 + 0.28209479177387814 * vec3<f32>(sogCode(sh0.r).g, sogCode(sh0.g).g, sogCode(sh0.b).g), 1.0);
        #endif
    #else
        let position = read_tPositionGroup(vertexId).xyz;
        // an integer (\`rgb32f-i32\`)
        let group = readGroup_tPositionGroup(vertexId);
    #endif

    ${assign_varyings_wgsl}

    #if defined(dSplatEncoding_sog)
        let q = normalize(sogRotation(vertexId));
        let sb = sogBytes(textureLoad(tSogScales, sogTexel(vertexId), 0)).rgb;
        let s = exp(vec3<f32>(sogCode(sb.r).r, sogCode(sb.g).r, sogCode(sb.b).r)) * uObject.uSizeFactor;
        #if defined(dSogVersion3)
            // to 1.5625 in the file: clamped (as \`gaussianSplatsFromSog\`)
            out.vSplatOpacity = min(1.0, (sh0.a / 204.0) * (sh0.a / 204.0));
        #else
            out.vSplatOpacity = sh0.a / 255.0;
        #endif
    #else
        // rotation (x, y, z snorm16, w unorm16) and log-scales (fixed point, 1/1024) with the opacity (unorm16)
        let r = read_tSplat(vertexId * 2);
        let so = read_tSplat(vertexId * 2 + 1);
        let q = normalize(vec4<f32>(r.xyz / 32767.0, r.w / 65535.0));
        let s = exp(so.xyz / 1024.0) * uObject.uSizeFactor;
        out.vSplatOpacity = so.w / 65535.0;
    #endif

    let xx = q.x * q.x;
    let yy = q.y * q.y;
    let zz = q.z * q.z;
    let xy = q.x * q.y;
    let xz = q.x * q.z;
    let yz = q.y * q.z;
    let wx = q.w * q.x;
    let wy = q.w * q.y;
    let wz = q.w * q.z;
    // R S: the rotation's columns scaled
    let rs = mat3x3<f32>(
        vec3<f32>(1.0 - 2.0 * (yy + zz), 2.0 * (xy + wz), 2.0 * (xz - wy)) * s.x,
        vec3<f32>(2.0 * (xy - wz), 1.0 - 2.0 * (xx + zz), 2.0 * (yz + wx)) * s.y,
        vec3<f32>(2.0 * (xz + wy), 2.0 * (yz - wx), 1.0 - 2.0 * (xx + yy)) * s.z
    );

    let aTransform = affineTransform(in.aTransform0, in.aTransform1, in.aTransform2);
    let world = uObject.uObjectTransform * (aTransform * vec4<f32>(position, 1.0));
    out.vModelPosition = (uGlobal.uModel * world).xyz;
    let modelView = uGlobal.uModelView * uObject.uObjectTransform * aTransform;
    let mvPosition = uGlobal.uModelView * world;
    out.vViewPosition = mvPosition.xyz;

    // from the camera (the view's origin) to the splat in its own frame: the model-view's rotation undone (its
    // transpose: rotations and uniform scales)
    let W = mat3x3<f32>(modelView[0].xyz, modelView[1].xyz, modelView[2].xyz);
    out.vSplatSh = vec3<f32>(0.0);
    #if defined(dSplatEncoding_sog) && defined(dSogSh)
        // the palette's entry by the splat's label (rg, version 3 rgb); a label past the palette none (as PlayCanvas)
        let lb = sogBytes(textureLoad(tSogShLabels, sogTexel(vertexId), 0));
        #if defined(dSogVersion3)
            let label = i32(lb.r + lb.g * 256.0 + lb.b * 65536.0);
        #else
            let label = i32(lb.r + lb.g * 256.0);
        #endif
        if (uObject.uSplatShDegree > 0 && label < uObject.uSogShCount) {
            out.vSplatSh = evalSplatSh(label, normalize(transpose(W) * mvPosition.xyz));
        }
    #elif !defined(dSplatEncoding_sog)
        if (uObject.uSplatShDegree > 0) {
            out.vSplatSh = evalSplatSh(vertexId, normalize(transpose(W) * mvPosition.xyz));
        }
    #endif

    let projection = uGlobal.uProjection;
    let viewport = uGlobal.uViewport;
    let ortho = projection[2][3] == 0.0;
    let focal = vec2<f32>(projection[0][0], projection[1][1]) * viewport.zw * 0.5;
    var J: mat3x3<f32>;
    if (ortho) {
        J = mat3x3<f32>(focal.x, 0.0, 0.0, 0.0, focal.y, 0.0, 0.0, 0.0, 0.0);
    } else {
        // the Jacobian of the perspective projection at the center (columns: d/dx, d/dy, d/dz of the pixel), the
        // center's direction clamped to 1.3 times the frustum (splats far off screen not stretched without bound)
        let tz = -mvPosition.z;
        let limit = 1.3 / vec2<f32>(projection[0][0], projection[1][1]);
        let t = clamp(mvPosition.xy / tz, -limit, limit) * tz;
        J = mat3x3<f32>(
            focal.x / tz, 0.0, 0.0,
            0.0, focal.y / tz, 0.0,
            focal.x * t.x / (tz * tz), focal.y * t.y / (tz * tz), 0.0
        );
    }
    let T = J * W * rs;
    let cov = T * transpose(T);
    // dilated by 0.3 px² against aliasing of splats smaller than a pixel
    let a = cov[0][0] + 0.3;
    let b = cov[0][1];
    let c = cov[1][1] + 0.3;

    // the 2D Gaussian's principal axes, in pixels, to uKernelSize standard deviations (at most 1024 px)
    let mid = 0.5 * (a + c);
    let radius = length(vec2<f32>(0.5 * (a - c), b));
    let lambda1 = mid + radius;
    let lambda2 = max(mid - radius, 0.1);
    var axis = select(vec2<f32>(0.0, 1.0), vec2<f32>(1.0, 0.0), a >= c);
    if (abs(b) > 1e-9) {
        axis = normalize(vec2<f32>(b, lambda1 - a));
    }
    let sigma = sqrt(vec2<f32>(lambda1, lambda2));
    let extent = min(uObject.uKernelSize * sigma, vec2<f32>(1024.0));
    let offset = mapping.x * extent.x * axis + mapping.y * extent.y * vec2<f32>(-axis.y, axis.x);
    out.vQuad = mapping * extent / sigma;

    let p = projection * mvPosition;
    out.position = vec4<f32>(p.xy + 2.0 * offset / viewport.zw * p.w, p.zw);

    if (!ortho && mvPosition.z >= 0.0) {
        // behind the camera: move out of [0, w] to 'discard' in the vertex stage
        out.position.z = 2.0 * out.position.w;
    }
    #if defined(dClipVariant_pixel) && defined(dClipObjects)
        // pixel clipping: a splat tested at its center; the ghost pass (\`uClipGhost\`) keeps what the clip objects remove
        if (clipTest(out.vModelPosition / uGlobal.uModelScale, clipFlag) != (uGlobal.uClipGhost != 0u)) {
            out.position.z = 2.0 * out.position.w;
        }
    #endif
    ${clip_instance_wgsl}
    return out;
}

@fragment
fn fragmentMain(in: VertexOutput) -> FragmentOutput {
    var out: FragmentOutput;

    if (fadeLod(in.position, in.vModelPosition)) {
        discard;
    }
    // pixel clipping in the vertex stage: a splat tested at its center

    let fragmentDepth = in.position.z;
    let interior = false;

    ${assign_material_wgsl}
    #if defined(dRenderVariant_color) || defined(dRenderVariant_tracing)
        // the view-dependent color over the base color, not below 0 (as 3D Gaussian splatting clamps)
        material = vec4<f32>(max(material.rgb + in.vSplatSh, vec3<f32>(0.0)), material.a);
    #endif

    // the Gaussian's falloff times the splat's opacity, in every variant (picked where it is seen)
    let coverage = exp(-0.5 * dot(in.vQuad, in.vQuad)) * in.vSplatOpacity;
    if (coverage < uObject.uSplatAlphaCutoff) {
        discard;
    }
    material.a *= coverage;
    #if defined(dRenderVariant_pick) || defined(dRenderVariant_marking)
        // picked and outlined where at least as opaque as the picking threshold: not by the faint tails of splats in
        // front of what is seen
        if (coverage < uGlobal.uPickingAlphaThreshold) {
            discard;
        }
    #endif

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
        var color = material;
        ${apply_marker_wgsl}
        #if defined(dRenderVariant_tracing)
            ${tracing_output_wgsl('normalize(in.vViewPosition)')}
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
