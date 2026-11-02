/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

/**
 * The geometry of a sphere on screen, free of any binding: shared by the sphere impostors (`spheres.wgsl.ts`) and
 * compute kernels that draw the same spheres
 */
export const sphere_shared_wgsl = /* wgsl */ `
/**
 * The bounding rectangle of a perspective-projected 3D sphere at view-space \`p\` with radius \`r\`: min x, min y, max x,
 * max y as \`sphereProjection\` combines them (scaled by \`projection\`'s focal terms, before its off-axis terms).
 * Michael Mara, Morgan McGuire. 2013; specialization by Arseny Kapoulkine, MIT License Copyright (c) 2018
 * https://github.com/zeux/niagara
 */
fn sphereBounds(p: vec3<f32>, r: f32, projection: mat4x4<f32>) -> vec4<f32> {
    let pr = p * r;
    let pzr2 = p.z * p.z - r * r;

    let vx = sqrt(p.x * p.x + pzr2);
    let minx = ((vx * p.x - pr.z) / (vx * p.z + pr.x)) * projection[0][0];
    let maxx = ((vx * p.x + pr.z) / (vx * p.z - pr.x)) * projection[0][0];

    let vy = sqrt(p.y * p.y + pzr2);
    let miny = ((vy * p.y - pr.z) / (vy * p.z + pr.y)) * projection[1][1];
    let maxy = ((vy * p.y + pr.z) / (vy * p.z - pr.y)) * projection[1][1];
    return vec4<f32>(minx, miny, maxx, maxy);
}

/**
 * The least and greatest \`x / z\` over the part beyond the near plane (\`z <= -near\`) of the circle at \`(a, c)\` of
 * radius \`r\` (a sphere's projection onto the plane of an axis and z): at the tangent points from the eye that lie beyond
 * it, and at the ends of the chord the plane cuts (the only bounds where the eye is in the circle); least above greatest
 * where no part lies beyond it. After "2D Polyhedral Bounds of a Clipped, Perspective-Projected 3D Sphere", Michael Mara,
 * Morgan McGuire. 2013
 */
fn nearAxisBounds(a: f32, c: f32, r: f32, near: f32) -> vec2<f32> {
    var lo = 1e30;
    var hi = -1e30;
    let l2 = a * a + c * c;
    let t2 = l2 - r * r;
    if (t2 > 0.0) {
        let t = sqrt(t2);
        for (var s = -1.0; s <= 1.0; s += 2.0) {
            let tx = (a * t2 - s * c * t * r) / l2;
            let tz = (c * t2 + s * a * t * r) / l2;
            if (tz <= -near) {
                lo = min(lo, tx / tz);
                hi = max(hi, tx / tz);
            }
        }
    }
    let d = c + near;
    let w2 = r * r - d * d;
    if (w2 > 0.0) {
        let w = sqrt(w2);
        lo = min(lo, min((a - w) / -near, (a + w) / -near));
        hi = max(hi, max((a - w) / -near, (a + w) / -near));
    }
    return vec2<f32>(lo, hi);
}

/**
 * \`sphereBounds\` of the part beyond the near plane of a sphere through it (for its bounds the full projection's are
 * unbounded where the sphere reaches behind the eye): min x above max x where no part lies beyond it
 */
fn sphereNearBounds(p: vec3<f32>, r: f32, near: f32, projection: mat4x4<f32>) -> vec4<f32> {
    let x = nearAxisBounds(p.x, p.z, r, near) * projection[0][0];
    let y = nearAxisBounds(p.y, p.z, r, near) * projection[1][1];
    return vec4<f32>(x.x, y.x, x.y, y.y);
}

/**
 * The pixels of a \`viewport\`-sized target the sphere of \`sphereBounds\` (perspective) spans: min x, min y, max x,
 * max y in pixel coordinates (y down), continuous
 */
fn spherePixelBounds(p: vec3<f32>, r: f32, projection: mat4x4<f32>, viewport: vec2<f32>) -> vec4<f32> {
    let b = sphereBounds(p, r, projection);
    // the bounds are negated, without the off-axis terms
    let x = (vec2<f32>(-b.z, -b.x) - projection[2].x) * 0.5 + 0.5;
    let y = (vec2<f32>(b.y, b.w) + projection[2].y) * 0.5 + 0.5;
    return vec4<f32>(x.x * viewport.x, y.x * viewport.y, x.y * viewport.x, y.y * viewport.y);
}

/**
 * Whether the sphere raster draws the sphere at view-space \`p\` with radius \`r\` (else the impostors do): in front of
 * the \`near\` plane as a whole and at most \`maxRadius\` pixels in half size on screen (its \`spherePixelBounds\`
 * \`bounds\`); both paths decide by it, so each sphere is drawn by one
 */
fn rasterizesSphere(p: vec3<f32>, r: f32, near: f32, bounds: vec4<f32>, maxRadius: f32) -> bool {
    let half = (bounds.zw - bounds.xy) * 0.5;
    return -p.z - r > near && max(half.x, half.y) <= maxRadius;
}

/**
 * A ray against a sphere of radius \`r\` whose center is at \`toCenter\` from the ray's origin, along unit \`direction\`:
 * x the discriminant (negative: no hit), y and z the ray parameters \`B - sqrt(det)\` and \`B + sqrt(det)\` (where it hits)
 */
fn raySphere(direction: vec3<f32>, toCenter: vec3<f32>, r: f32) -> vec3<f32> {
    let B = dot(direction, toCenter);
    let det = B * B + r * r - dot(toCenter, toCenter);
    let s = sqrt(max(det, 0.0));
    return vec3<f32>(det, B - s, B + s);
}
`;
