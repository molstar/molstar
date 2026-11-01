/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

/** The signed distance of a point in a clip object's frame to its shape (\`localSD\`) */
export const clip_shape_wgsl = /* wgsl */ `
/**
 * Signed distance of \`t\` (in the shape's frame) to a shape of half size \`size\` (1 plane, 2 sphere, 3 cube,
 * 4 cylinder, 5 infinite cone)
 */
fn localSD(type_: i32, t: vec3<f32>, size: vec3<f32>) -> f32 {
    if (type_ == 1) {
        return -t.y;
    } else if (type_ == 2) {
        return (length(t / size) - 1.0) * min(min(size.x, size.y), size.z);
    } else if (type_ == 3) {
        let d = abs(t) - size;
        return min(max(d.x, max(d.y, d.z)), 0.0) + length(max(d, vec3<f32>(0.0)));
    } else if (type_ == 4) {
        let d = abs(vec2<f32>(length(t.xz), t.y)) - size.xy;
        return min(max(d.x, d.y), 0.0) + length(max(d, vec2<f32>(0.0)));
    } else if (type_ == 5) {
        return dot(size.xy, vec2<f32>(length(t.xy), t.z));
    }
    return 0.1;
}
`;

/**
 * What clips by clip objects (`localSD` of `clip_shape_wgsl`): the signed distance, `clipTest`, `clipDistance`,
 * `clipNormal`, the intervals along a segment, `clipExit` and a sphere's class (`clipSphereClass`), on accessors of clip object `i` its caller defines
 * (`clipObjectFrame`, `clipObjectShape`, `isClipObjectActive(i, flag)`): the programs' from the material block, the
 * kernels' from storage (`clip_storage_wgsl`); `count` the clip objects (the programs' a pipeline constant, the
 * loops unrolled)
 */
function clip_core_wgsl(count: string) {
  return /* wgsl */ `
/** The point \`p\` in the frame of clip object \`i\` (\`Clip.Objects.frame\`) */
fn clipLocal(i: i32, p: vec3<f32>) -> vec3<f32> {
    return (clipObjectFrame(i) * vec4<f32>(p, 1.0)).xyz;
}

/** The type of clip object \`i\`: \`w\` of its shape, negative where inverted (\`Clip.Objects.shape\`) */
fn clipType(i: i32) -> i32 {
    return i32(abs(clipObjectShape(i).w) + 0.5);
}

fn getSignedDistance(center: vec3<f32>, i: i32) -> f32 {
    return localSD(clipType(i), clipLocal(i, center), clipObjectShape(i).xyz);
}

fn isClipObjectInverted(i: i32) -> bool {
    return clipObjectShape(i).w < 0.0;
}

fn clipTest(center: vec3<f32>, flag: i32) -> bool {
    for (var i = 0; i < ${count}; i++) {
        if (isClipObjectActive(i, flag)) {
            let test = getSignedDistance(center, i) <= 0.0;
            if (test != isClipObjectInverted(i)) {
                return true;
            }
        }
    }
    return false;
}

fn clipDistance(center: vec3<f32>, flag: i32) -> f32 {
    var d = 1.0e10;
    for (var i = 0; i < ${count}; i++) {
        if (isClipObjectActive(i, flag)) {
            let sd = getSignedDistance(center, i);
            d = min(d, select(sd, -sd, isClipObjectInverted(i)));
        }
    }
    return d;
}

fn clipNormal(center: vec3<f32>, flag: i32) -> vec3<f32> {
    let e = vec2<f32>(0.01, 0.0);
    return normalize(vec3<f32>(
        clipDistance(center + e.xyy, flag) - clipDistance(center - e.xyy, flag),
        clipDistance(center + e.yxy, flag) - clipDistance(center - e.yxy, flag),
        clipDistance(center + e.yyx, flag) - clipDistance(center - e.yyx, flag)
    ));
}

fn clipHalfLine(f0: f32, f1: f32) -> vec2<f32> {
    let df = f1 - f0;
    if (abs(df) < 1.0e-12) {
        return select(vec2<f32>(1.0e6, -1.0e6), vec2<f32>(-1.0e6, 1.0e6), f0 <= 0.0);
    }
    let sr = -f0 / df;
    return select(vec2<f32>(-1.0e6, sr), vec2<f32>(sr, 1.0e6), df < 0.0);
}

fn clipSlab(q0: f32, q1: f32, size: f32) -> vec2<f32> {
    let a = clipHalfLine(q0 - size, q1 - size);
    let b = clipHalfLine(-q0 - size, -q1 - size);
    return vec2<f32>(max(a.x, b.x), min(a.y, b.y));
}

fn clipQuadratic(a: f32, b: f32, c: f32) -> vec2<f32> {
    if (a < 1.0e-12) {
        return select(vec2<f32>(1.0e6, -1.0e6), vec2<f32>(-1.0e6, 1.0e6), c <= 0.0);
    }
    let disc = b * b - a * c;
    if (disc < 0.0) {
        return vec2<f32>(1.0e6, -1.0e6);
    }
    let sq = sqrt(disc);
    return vec2<f32>((-b - sq) / a, (-b + sq) / a);
}

/** Interval of the segment from \`pFront\` (0) to \`pBack\` (1) inside clip object \`i\` */
fn clipObjectInterval(i: i32, pFront: vec3<f32>, pBack: vec3<f32>) -> vec2<f32> {
    let type_ = clipType(i);
    let qF = clipLocal(i, pFront);
    let qB = clipLocal(i, pBack);
    if (type_ == 1) {
        return clipHalfLine(-qF.y, -qB.y);
    }
    let size = clipObjectShape(i).xyz;
    if (type_ == 2) {
        let o = qF / size;
        let d = (qB - qF) / size;
        return clipQuadratic(dot(d, d), dot(o, d), dot(o, o) - 1.0);
    } else if (type_ == 3) {
        let x = clipSlab(qF.x, qB.x, size.x);
        let y = clipSlab(qF.y, qB.y, size.y);
        let z = clipSlab(qF.z, qB.z, size.z);
        return vec2<f32>(max(x.x, max(y.x, z.x)), min(x.y, min(y.y, z.y)));
    } else if (type_ == 4) {
        let o = qF.xz;
        let d = (qB - qF).xz;
        let body = clipQuadratic(dot(d, d), dot(o, d), dot(o, o) - size.x * size.x);
        let slab = clipSlab(qF.y, qB.y, size.y);
        return vec2<f32>(max(body.x, slab.x), min(body.y, slab.y));
    } else if (type_ == 5) {
        let xyF = qF.xy;
        let xyd = (qB - qF).xy;
        let hF = -size.y * qF.z;
        let hd = -size.y * (qB.z - qF.z);
        let nappe = clipHalfLine(-hF, -hF - hd);
        let a = size.x * size.x * dot(xyd, xyd) - hd * hd;
        let b = size.x * size.x * dot(xyF, xyd) - hF * hd;
        let c = size.x * size.x * dot(xyF, xyF) - hF * hF;
        if (a > 1.0e-12) {
            let g = clipQuadratic(a, b, c);
            return vec2<f32>(max(g.x, nappe.x), min(g.y, nappe.y));
        }
        let disc = b * b - a * c;
        if (a < -1.0e-12 && disc >= 0.0) {
            let sq = sqrt(disc);
            let s1 = (-b + sq) / a;
            let s2 = (-b - sq) / a;
            let lo = vec2<f32>(max(-1.0e6, nappe.x), min(s1, nappe.y));
            let hi = vec2<f32>(max(s2, nappe.x), min(1.0e6, nappe.y));
            return select(hi, lo, lo.x <= lo.y);
        }
        let lin = clipHalfLine(c, 2.0 * b + a + c);
        return vec2<f32>(max(lin.x, nappe.x), min(lin.y, nappe.y));
    }
    return vec2<f32>(1.0e6, -1.0e6);
}

/**
 * First position from \`sStart\` along the segment from \`pFront\` (0) to \`pBack\` (1) that no clip object removes,
 * -1 if there is none
 */
fn clipExit(pFront: vec3<f32>, pBack: vec3<f32>, sStart: f32, flag: i32) -> f32 {
    var s = sStart;
    for (var k = 0; k <= ${count}; k++) {
        var sNext = s;
        for (var i = 0; i < ${count}; i++) {
            if (isClipObjectActive(i, flag)) {
                let iv = clipObjectInterval(i, pFront, pBack);
                if (isClipObjectInverted(i)) {
                    if (iv.x > iv.y || s > iv.y) {
                        sNext = 1.0e6;
                    } else if (s < iv.x) {
                        sNext = max(sNext, iv.x);
                    }
                } else if (s > iv.x && s < iv.y) {
                    sNext = max(sNext, iv.y);
                }
            }
        }
        if (sNext == s) {
            return s;
        }
        s = sNext;
        if (s >= 1.0) {
            return -1.0;
        }
    }
    return -1.0;
}

/**
 * The sphere at \`p\` (the clip objects' space) of radius \`r\` against the clip objects, from the signed distance of its center: 0 no
 * point clipped, 2 every point (by one object), 1 per pixel. The distances (\`localSD\`) are at most the true ones once
 * divided by their largest slope (1, the cone's \`|size.xy|\`; the ellipsoid's scaled by its smallest axis): bounds, the
 * radius by the frame's largest scale
 */
fn clipSphereClass(p: vec3<f32>, r: f32, flag: i32) -> u32 {
    var perPixel = false;
    for (var i = 0; i < ${count}; i++) {
        if (!isClipObjectActive(i, flag) || clipType(i) == 0) { continue; }
        let f = clipObjectFrame(i);
        var rl = r * max(length(f[0].xyz), max(length(f[1].xyz), length(f[2].xyz)));
        if (clipType(i) == 5) {
            rl *= length(clipObjectShape(i).xy);
        }
        let sd = getSignedDistance(p, i);
        // within the clipped side at least this far (it is \`sd <= 0\`, inverted \`sd > 0\`)
        let depth = select(-sd, sd, isClipObjectInverted(i));
        if (depth > rl) {
            return 2u;
        }
        if (depth >= -rl) {
            perPixel = true;
        }
    }
    return select(0u, 1u, perPixel);
}
`;
}

/**
 * Clipping by clip objects read from storage (compute kernels, which have no material block): per clip object five
 * vec4s in `buffer` from vec4 `base`, its frame's columns (`Clip.Objects.frame`) and its shape (half size, the type
 * in `w`, negative where inverted), `count` of them (expressions); the functions of `clip_core_wgsl`, without a flag
 * (no clipping layer: pass 0) or a solid interior's clip object
 */
export function clip_storage_wgsl(buffer: string, base: string, count: string) {
  return /* wgsl */ `
${clip_shape_wgsl}

fn clipObjectFrame(i: i32) -> mat4x4<f32> {
    let o = ${base} + u32(i) * 5u;
    return mat4x4<f32>(${buffer}[o], ${buffer}[o + 1u], ${buffer}[o + 2u], ${buffer}[o + 3u]);
}

fn clipObjectShape(i: i32) -> vec4<f32> {
    return ${buffer}[${base} + u32(i) * 5u + 4u];
}

fn isClipObjectActive(i: i32, flag: i32) -> bool {
    return flag == 0 || (flag & (1 << u32(i))) == 0;
}

${clip_core_wgsl(`i32(${count})`)}
`;
}

/**
 * Clip object signed distances, tests and ray intervals, bounded by \`oClipObjectCount\`. \`flag\` is the bit-flag of the
 * clip objects to ignore (object ids start at 1), see \`clip_flag_wgsl\`.
 */
export const clip_functions_wgsl = /* wgsl */ `
fn quaternionTransform(q: vec4<f32>, v: vec3<f32>) -> vec3<f32> {
    let t = 2.0 * cross(q.xyz, v);
    return v + q.w * t + cross(q.xyz, t);
}

${clip_shape_wgsl}

/** Signed distance of \`center\` to a shape by position, rotation, scale and the transform of the point (image trims) */
fn getShapeSignedDistance(center: vec3<f32>, type_: i32, position: vec3<f32>, rotation: vec4<f32>, scale: vec3<f32>, transform: mat4x4<f32>) -> f32 {
    let c = (transform * vec4<f32>(center, 1.0)).xyz;
    return localSD(type_, quaternionTransform(vec4<f32>(-rotation.xyz, rotation.w), c - position), scale * 0.5);
}

#if defined(dClipObjects)
override oClipObjectCount: i32 = 1;

/** The clip objects of the program (\`uMaterial\`), the core's accessors */
fn clipObjectFrame(i: i32) -> mat4x4<f32> {
    return uMaterial.uClipObjectFrame[i];
}

fn clipObjectShape(i: i32) -> vec4<f32> {
    return uMaterial.uClipObjectShape[i];
}

fn isClipObjectActive(i: i32, flag: i32) -> bool {
    return i != uDraw.uSolidInteriorClip && (flag == 0 || (flag & (1 << u32(i))) == 0);
}

${clip_core_wgsl('oClipObjectCount')}

/** Whether the clip object \`index\` doesn't cut the fragment (its group is exempt) */
fn clipIgnored(index: i32, flag: i32) -> bool {
    return flag != 0 && (flag & (1 << u32(index))) != 0;
}

fn isSolidInteriorClip(i: i32) -> bool {
    return i >= 0 && i < oClipObjectCount;
}

/**
 * Where the segment from \`pFront\` (0) to \`pBack\` (1) enters the solid interior clip object (\`uSolidInteriorClip\`),
 * or leaves it for \`back\`, -1 if it misses it
 */
fn clipCapExit(pFront: vec3<f32>, pBack: vec3<f32>, back: bool) -> f32 {
    let i = uDraw.uSolidInteriorClip;
    if (!isSolidInteriorClip(i)) {
        return -1.0;
    }
    let iv = clipObjectInterval(i, pFront, pBack);
    if (iv.x < iv.y) {
        return select(iv.y, iv.x, isClipObjectInverted(i) != back);
    }
    return -1.0;
}

fn clipCapDistance(center: vec3<f32>) -> f32 {
    let i = uDraw.uSolidInteriorClip;
    if (!isSolidInteriorClip(i)) {
        return 0.0;
    }
    let sd = getSignedDistance(center, i);
    return select(sd, -sd, isClipObjectInverted(i));
}

fn clipCapNormal(center: vec3<f32>) -> vec3<f32> {
    let e = vec2<f32>(0.01, 0.0);
    return normalize(vec3<f32>(
        clipCapDistance(center + e.xyy) - clipCapDistance(center - e.xyy),
        clipCapDistance(center + e.yxy) - clipCapDistance(center - e.yxy),
        clipCapDistance(center + e.yyx) - clipCapDistance(center - e.yyx)
    ));
}
#endif
`;

/** The clip object bit-flag \`clipFlag\` from \`clipping\` (an expression of the clipping varying), 0 without one */
export function clip_flag_wgsl(clipping: string) {
  return /* wgsl */ `
#if defined(dClipObjects) && defined(dClipping)
    let clipFlag = i32(floor(${clipping} * 255.0 + 0.5));
#else
    let clipFlag = 0;
#endif
`;
}

/** Instance clipping in the vertex stage: moves the vertex out of the clip volume; expects \`aTransform\` and \`clipFlag\` */
export const clip_instance_wgsl = /* wgsl */ `
#if defined(dClipVariant_instance) && defined(dClipObjects)
    let clipCenter = (uGlobal.uModel * aTransform * vec4<f32>(invariantBoundingSphere(in).xyz, 1.0)).xyz;
    if (clipTest(clipCenter / uGlobal.uModelScale, clipFlag)) {
        // move out of [0, w] to 'discard' in the vertex stage
        out.position.z = 2.0 * out.position.w;
    }
#endif
`;

/** Pixel clipping of \`modelPosition\` (an expression) in the fragment stage; expects \`clipFlag\` */
export function clip_pixel_wgsl(modelPosition: string) {
  return /* wgsl */ `
#if defined(dClipVariant_pixel) && defined(dClipObjects)
    if (clipTest(${modelPosition} / uGlobal.uModelScale, clipFlag)) {
        discard;
    }
#endif
`;
}
