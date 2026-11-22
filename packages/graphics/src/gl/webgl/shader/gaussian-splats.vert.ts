/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { common } from './chunks/common.glsl.js';
import { read_from_texture } from './chunks/read-from-texture.glsl.js';
import { common_vert_params } from './chunks/common-vert-params.glsl.js';
import { color_vert_params } from './chunks/color-vert-params.glsl.js';
import { common_clip } from './chunks/common-clip.glsl.js';
import { assign_color_varying } from './chunks/assign-color-varying.glsl.js';
import { assign_marker_varying } from './chunks/assign-marker-varying.glsl.js';
import { assign_clipping_varying } from './chunks/assign-clipping-varying.glsl.js';
import { clip_instance } from './chunks/clip-instance.glsl.js';

/**
 * Gaussian splats as quads of two triangles (6 vertices) along the principal axes of each splat's projected 2D
 * Gaussian (EWA): the 3D covariance `R S Sᵀ Rᵀ` from the splat's rotation and scales (`tSplat`), projected by the
 * Jacobian of the projection at its center
 */
export const gaussianSplats_vert = `
precision highp float;
precision highp int;

${common}
${read_from_texture}
${common_vert_params}
${color_vert_params}
${common_clip}

/** the fragment's offset from the center in standard deviations along the 2D Gaussian's axes */
varying vec2 vQuad;
varying float vSplatOpacity;
/** the view-dependent color past the base color (\`tSplatSh\`), zero without */
varying vec3 vSplatSh;

#if defined(dSplatEncoding_sog)
    /** the bytes of texel \`p\` of SOG's image \`t\` of \`dim\` texels */
    vec4 sogTexel(const in sampler2D t, const in vec2 p, const in vec2 dim) {
        #if __VERSION__ == 100
            vec4 v = texture2D(t, (p + 0.5) / dim);
        #else
            vec4 v = texelFetch(t, ivec2(p), 0);
        #endif
        return floor(v * 255.0 + 0.5);
    }

    /** the bytes of splat \`i\` of SOG's image \`t\` (a pixel a splat, row by row) */
    vec4 sogSplat(const in sampler2D t, const in float i) {
        float y = floor((i + 0.5) / uSogImageDim.x);
        return sogTexel(t, vec2(i - y * uSogImageDim.x, y), uSogImageDim);
    }

    /** the codebooks' entry \`b\`: scale (r), base color (g), view-dependent color (b) */
    vec4 sogCode(const in float b) {
        #if __VERSION__ == 100
            return texture2D(tSogCodebook, vec2((b + 0.5) / 256.0, 0.5));
        #else
            return texelFetch(tSogCodebook, ivec2(int(b), 0), 0);
        #endif
    }

    /** splat \`i\`'s center: 16 bits an axis within the means' range, out of the log domain */
    vec3 sogPosition(const in float i) {
        vec3 lo = sogSplat(tSogMeansLow, i).rgb;
        vec3 hi = sogSplat(tSogMeansHigh, i).rgb;
        vec3 n = uSogMeansMin + (uSogMeansMax - uSogMeansMin) * ((lo + hi * 256.0) / 65535.0);
        return sign(n) * (exp(abs(n)) - 1.0);
    }

    /** splat \`i\`'s rotation (\`x, y, z, w\`): version 3 an exponential map, version 2 the smallest three */
    vec4 sogRotation(const in float i) {
        vec4 b = sogSplat(tSogQuats, i);
        #if defined(dSogVersion3)
            vec3 p = (b.xyz - 128.0) / 127.0;
            float r = length(p);
            float t = 1.5707963267948966 * r;
            float k = r > 1e-12 ? sin(t) / r : 1.5707963267948966;
            return vec4(p * k, cos(t));
        #else
            // the largest component's index in \`w, x, y, z\` (252 to 255), the stored three the others in order
            float largest = b.w - 252.0;
            vec3 c = (b.xyz / 255.0 * 2.0 - 1.0) * 0.7071067811865476;
            float m = sqrt(max(0.0, 1.0 - dot(c, c)));
            if (largest < -0.5 || largest > 3.5) return vec4(0.0, 0.0, 0.0, 1.0);
            if (largest < 0.5) return vec4(c, m);
            if (largest < 1.5) return vec4(m, c.y, c.z, c.x);
            if (largest < 2.5) return vec4(c.y, m, c.z, c.x);
            return vec4(c.y, c.z, m, c.x);
        #endif
    }
#endif

#if defined(dSplatEncoding_sog)
    /**
     * the palette's coefficient \`k\` of entry \`label\` (rgb, \`uSogShCoefficients\` pixels an entry, 64 entries a row),
     * see \`getSogShLayout\`; \`count\` unused (the palette's own)
     */
    vec3 splatSh(const in int label, const in int count, const in int k) {
        #if defined(dSogSh)
            float u = mod(float(label), 64.0) * float(uSogShCoefficients) + float(k);
            vec3 b = sogTexel(tSogShCentroids, vec2(u, floor(float(label) / 64.0)), uSogShDim).rgb;
            return vec3(sogCode(b.r).b, sogCode(b.g).b, sogCode(b.b).b);
        #else
            return vec3(0.0);
        #endif
    }
#else
    /** splat \`splat\`'s coefficient of basis function \`k\` past the zeroth (rgb), see \`SplatSh\` */
    vec3 splatSh(const in int splat, const in int count, const in int k) {
        return (read_tSplatSh(splat * count + k).rgb * 2.0 - 1.0) * uSplatShScale;
    }
#endif

/**
 * The view-dependent color of splat \`splat\` past its base color (3D Gaussian splatting's \`eval_sh\` without the zeroth),
 * \`dir\` the unit direction from the camera to it in its own frame
 */
vec3 evalSplatSh(const in int splat, const in vec3 dir) {
    int degree = uSplatShDegree;
    int n = degree == 1 ? 3 : (degree == 2 ? 8 : 15);
    float x = dir.x, y = dir.y, z = dir.z;
    vec3 c = 0.4886025119029199 * (-y * splatSh(splat, n, 0) + z * splatSh(splat, n, 1) - x * splatSh(splat, n, 2));
    if (degree > 1) {
        float xx = x * x, yy = y * y, zz = z * z;
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

void main(){
    #if defined(dInstanceData)
        readInstanceData();
    #endif
    vec2 mapping = vec2(1.0, 1.0); // vertices 2 and 5
    #if __VERSION__ == 100
        int m = imod(VertexID, 6);
    #else
        int m = VertexID % 6;
    #endif
    if (m == 0) {
        mapping = vec2(-1.0, 1.0);
    } else if (m == 1 || m == 3) {
        mapping = vec2(-1.0, -1.0);
    } else if (m == 4) {
        mapping = vec2(1.0, -1.0);
    }
    int vertexId = VertexID / 6;
    // the splat of this draw slot: the instance's block of the order, back to front for the view (\`tSplatOrder\`, a
    // u32 in four bytes)
    vec4 orderBytes = read_tSplatOrder(int(aInstance) * uVertexCount + vertexId) * 255.0 + 0.5;
    vertexId = int(orderBytes.x) + int(orderBytes.y) * 256 + int(orderBytes.z) * 65536 + int(orderBytes.w) * 16777216;
    #if defined(dSplatEncoding_sog)
        // SOG's images read in place (\`GaussianSplats.Sog\`): the group the splat's index
        float splat = float(vertexId);
        vec3 position = sogPosition(splat);
        GroupId group = GroupId(vertexId);
        vec4 sh0 = sogSplat(tSogSh0, splat);
        #if defined(dColorType_direct) && (defined(dRenderVariant_color) || defined(dRenderVariant_tracing) || defined(dRenderVariant_emissive))
            // the base color the file's (\`0.5 + SH_C0 dc\`), before overpaint mixes over it
            vColor = vec4(0.5 + 0.28209479177387814 * vec3(sogCode(sh0.r).g, sogCode(sh0.g).g, sogCode(sh0.b).g), 1.0);
        #endif
    #else
        vec3 position = read_tPositionGroup(vertexId).rgb;
        GroupId group = readGroup_tPositionGroup(vertexId);
    #endif

    ${assign_color_varying}
    ${assign_marker_varying}
    ${assign_clipping_varying}

    #if defined(dSplatEncoding_sog)
        vec4 q = normalize(sogRotation(splat));
        vec3 sb = sogSplat(tSogScales, splat).rgb;
        vec3 s = exp(vec3(sogCode(sb.r).r, sogCode(sb.g).r, sogCode(sb.b).r)) * uSizeFactor;
        #if defined(dSogVersion3)
            // to 1.5625 in the file: clamped (as \`gaussianSplatsFromSog\`)
            vSplatOpacity = min(1.0, (sh0.a / 204.0) * (sh0.a / 204.0));
        #else
            vSplatOpacity = sh0.a / 255.0;
        #endif
    #else
        // rotation (x, y, z snorm16, w unorm16) and log-scales (fixed point, 1/1024) with the opacity (unorm16)
        vec4 r = read_tSplat(vertexId * 2);
        vec4 so = read_tSplat(vertexId * 2 + 1);
        vec4 q = normalize(vec4(r.xyz / 32767.0, r.w / 65535.0));
        vec3 s = exp(so.xyz / 1024.0) * uSizeFactor;
        vSplatOpacity = so.w / 65535.0;
    #endif

    float xx = q.x * q.x, yy = q.y * q.y, zz = q.z * q.z;
    float xy = q.x * q.y, xz = q.x * q.z, yz = q.y * q.z;
    float wx = q.w * q.x, wy = q.w * q.y, wz = q.w * q.z;
    // R S: the rotation's columns scaled
    mat3 rs = mat3(
        vec3(1.0 - 2.0 * (yy + zz), 2.0 * (xy + wz), 2.0 * (xz - wy)) * s.x,
        vec3(2.0 * (xy - wz), 1.0 - 2.0 * (xx + zz), 2.0 * (yz + wx)) * s.y,
        vec3(2.0 * (xz + wy), 2.0 * (yz - wx), 1.0 - 2.0 * (xx + yy)) * s.z
    );

    vec4 world = uObjectTransform * (aTransform * vec4(position, 1.0));
    vModelPosition = (uModel * world).xyz; // for clipping in frag shader
    mat4 modelView = uModelView * uObjectTransform * aTransform;
    vec4 mvPosition = uModelView * world;
    vViewPosition = mvPosition.xyz;

    // from the camera (the view's origin) to the splat in its own frame: the model-view's rotation undone (its
    // transpose: rotations and uniform scales)
    mat3 W = mat3(modelView[0].xyz, modelView[1].xyz, modelView[2].xyz);
    #if defined(dSplatEncoding_sog) && defined(dSogSh)
        // the palette's entry by the splat's label (rg, version 3 rgb); a label past the palette none (as PlayCanvas)
        vec4 lb = sogSplat(tSogShLabels, splat);
        #if defined(dSogVersion3)
            int label = int(lb.r + lb.g * 256.0 + lb.b * 65536.0);
        #else
            int label = int(lb.r + lb.g * 256.0);
        #endif
        vSplatSh = uSplatShDegree > 0 && label < uSogShCount
            ? evalSplatSh(label, normalize(transpose3(W) * mvPosition.xyz))
            : vec3(0.0);
    #elif defined(dSplatEncoding_sog)
        vSplatSh = vec3(0.0);
    #else
        vSplatSh = uSplatShDegree > 0 ? evalSplatSh(vertexId, normalize(transpose3(W) * mvPosition.xyz)) : vec3(0.0);
    #endif

    bool ortho = uProjection[2][3] == 0.0;
    vec2 focal = vec2(uProjection[0][0], uProjection[1][1]) * uViewport.zw * 0.5;
    mat3 J;
    if (ortho) {
        J = mat3(focal.x, 0.0, 0.0, 0.0, focal.y, 0.0, 0.0, 0.0, 0.0);
    } else {
        // the Jacobian of the perspective projection at the center (columns: d/dx, d/dy, d/dz of the pixel), the
        // center's direction clamped to 1.3 times the frustum (splats far off screen not stretched without bound)
        float tz = -mvPosition.z;
        vec2 limit = 1.3 / vec2(uProjection[0][0], uProjection[1][1]);
        vec2 t = clamp(mvPosition.xy / tz, -limit, limit) * tz;
        J = mat3(
            focal.x / tz, 0.0, 0.0,
            0.0, focal.y / tz, 0.0,
            focal.x * t.x / (tz * tz), focal.y * t.y / (tz * tz), 0.0
        );
    }
    mat3 T = J * W * rs;
    mat3 cov = T * transpose3(T);
    // dilated by 0.3 px² against aliasing of splats smaller than a pixel
    float a = cov[0][0] + 0.3;
    float b = cov[0][1];
    float c = cov[1][1] + 0.3;

    // the 2D Gaussian's principal axes, in pixels, to uKernelSize standard deviations (at most 1024 px)
    float mid = 0.5 * (a + c);
    float radius = length(vec2(0.5 * (a - c), b));
    float lambda1 = mid + radius;
    float lambda2 = max(mid - radius, 0.1);
    vec2 axis = abs(b) > 1e-9 ? normalize(vec2(b, lambda1 - a)) : (a >= c ? vec2(1.0, 0.0) : vec2(0.0, 1.0));
    vec2 sigma = sqrt(vec2(lambda1, lambda2));
    vec2 extent = min(uKernelSize * sigma, vec2(1024.0));
    vec2 offset = mapping.x * extent.x * axis + mapping.y * extent.y * vec2(-axis.y, axis.x);
    vQuad = mapping * extent / sigma;

    gl_Position = uProjection * mvPosition;
    gl_Position.xy += 2.0 * offset / uViewport.zw * gl_Position.w;

    if (!ortho && mvPosition.z >= 0.0) {
        // behind the camera: move out of [ -w, +w ] to 'discard' in vert shader
        gl_Position.z = 2.0 * gl_Position.w;
    }

    #if defined(dClipVariant_pixel) && defined(dClipObjects)
        // pixel clipping: a splat tested at its center; the ghost pass (\`uClipGhost\`) keeps what the clip objects remove
        if (clipTest(vModelPosition / uModelScale) != uClipGhost) {
            gl_Position.z = 2.0 * gl_Position.w;
        }
    #endif

    ${clip_instance}
}
`;
