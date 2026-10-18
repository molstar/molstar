/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

/**
 * The WGSL counterpart of `shader.frag.ts`, for grid3d compute kernels: the one-row textures (centers, info, coefficients,
 * alpha) are read by index with `textureLoad`, as the GLSL reads them at `index / (count - 1)` with nearest filtering
 */
export const UTILS_WGSL = /* wgsl */ `
fn L1(p: vec3<f32>, a0: f32, a1: f32, a2: f32) -> f32 {
    return a0 * p.z + a1 * p.x + a2 * p.y;
}

fn L2(p: vec3<f32>, a0: f32, a1: f32, a2: f32, a3: f32, a4: f32) -> f32 {
    let x = p.x;
    let y = p.y;
    let z = p.z;
    let xx = x * x;
    let yy = y * y;
    let zz = z * z;
    return (
        a0 * (-0.5 * xx - 0.5 * yy + zz) +
        a1 * (1.7320508075688772 * x * z) +
        a2 * (1.7320508075688772 * y * z) +
        a3 * (0.8660254037844386 * xx - 0.8660254037844386 * yy) +
        a4 * (1.7320508075688772 * x * y)
    );
}

fn L3(p: vec3<f32>, a0: f32, a1: f32, a2: f32, a3: f32, a4: f32, a5: f32, a6: f32) -> f32 {
    let x = p.x;
    let y = p.y;
    let z = p.z;
    let xx = x * x;
    let yy = y * y;
    let zz = z * z;
    let xxx = xx * x;
    let yyy = yy * y;
    let zzz = zz * z;
    return (
        a0 * (-1.5 * xx * z - 1.5 * yy * z + zzz) +
        a1 * (-0.6123724356957945 * xxx - 0.6123724356957945 * x * yy + 2.449489742783178 * x * zz) +
        a2 * (-0.6123724356957945 * xx * y - 0.6123724356957945 * yyy + 2.449489742783178 * y * zz) +
        a3 * (1.9364916731037085 * xx * z - 1.9364916731037085 * yy * z) +
        a4 * (3.872983346207417 * x * y * z) +
        a5 * (0.7905694150420949 * xxx - 2.3717082451262845 * x * yy) +
        a6 * (2.3717082451262845 * xx * y - 0.7905694150420949 * yyy)
    );
}

fn L4(p: vec3<f32>, a0: f32, a1: f32, a2: f32, a3: f32, a4: f32, a5: f32, a6: f32, a7: f32, a8: f32) -> f32 {
    let x = p.x;
    let y = p.y;
    let z = p.z;
    let xx = x * x;
    let yy = y * y;
    let zz = z * z;
    let xxx = xx * x;
    let yyy = yy * y;
    let zzz = zz * z;
    let xxxx = xxx * x;
    let yyyy = yyy * y;
    let zzzz = zzz * z;
    return (
        a0 * (0.375 * xxxx + 0.75 * xx * yy + 0.375 * yyyy - 3.0 * xx * zz - 3.0 * yy * zz + zzzz) +
        a1 * (-2.3717082451262845 * xxx * z - 2.3717082451262845 * x * yy * z + 3.1622776601683795 * x * zzz) +
        a2 * (-2.3717082451262845 * xx * y * z - 2.3717082451262845 * yyy * z + 3.1622776601683795 * y * zzz) +
        a3 * (-0.5590169943749475 * xxxx + 0.5590169943749475 * yyyy + 3.3541019662496847 * xx * zz - 3.3541019662496847 * yy * zz) +
        a4 * (-1.118033988749895 * xxx * y - 1.118033988749895 * x * yyy + 6.708203932499369 * x * y * zz) +
        a5 * (2.091650066335189 * xxx * z + -6.274950199005566 * x * yy * z) +
        a6 * (6.274950199005566 * xx * y * z + -2.091650066335189 * yyy * z) +
        a7 * (0.739509972887452 * xxxx - 4.437059837324712 * xx * yy + 0.739509972887452 * yyyy) +
        a8 * (2.958039891549808 * xxx * y + -2.958039891549808 * x * yyy)
    );
}

fn alpha(offset: i32) -> f32 {
    return textureLoad(tAlpha, vec2<i32>(offset, 0), 0).r;
}

fn Y(L: i32, X: vec3<f32>, aO: i32) -> f32 {
    if (L == 0) {
        return alpha(aO);
    } else if (L == 1) {
        return L1(X, alpha(aO), alpha(aO + 1), alpha(aO + 2));
    } else if (L == 2) {
        return L2(X, alpha(aO), alpha(aO + 1), alpha(aO + 2), alpha(aO + 3), alpha(aO + 4));
    } else if (L == 3) {
        return L3(X, alpha(aO), alpha(aO + 1), alpha(aO + 2), alpha(aO + 3), alpha(aO + 4), alpha(aO + 5), alpha(aO + 6));
    } else if (L == 4) {
        return L4(
            X,
            alpha(aO), alpha(aO + 1), alpha(aO + 2), alpha(aO + 3), alpha(aO + 4),
            alpha(aO + 5), alpha(aO + 6), alpha(aO + 7), alpha(aO + 8)
        );
    }
    return 0.0;
}

fn R(R2: f32, start: i32, end: i32) -> f32 {
    var gauss = 0.0;
    for (var i = start; i < end; i++) {
        let c = textureLoad(tCoeff, vec2<i32>(i, 0), 0).xy;
        gauss += c.x * exp(-c.y * R2);
    }
    return gauss;
}
`;

export const MAIN_WGSL = /* wgsl */ `
    var v = 0.0;
    for (var i = 0; i < uParams.uNCenters; i++) {
        let center = textureLoad(tCenters, vec2<i32>(i, 0), 0);
        let X = xyz - center.xyz;
        let R2 = dot(X, X);

        // center.w is the squared cutoff radius
        if (R2 > center.w) {
            continue;
        }

        let info = textureLoad(tInfo, vec2<i32>(i, 0), 0);
        v += R(R2, i32(info.z), i32(info.w)) * Y(i32(info.x), X, i32(info.y));
    }
`;
