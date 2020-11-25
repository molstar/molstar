/**
 * Copyright (c) 2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Inspired by https://github.com/dgasmith/gau2grid.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// gaussian:
//     R_0, R^+_1, R^-_1, ..., R^+_l, R^-_l
// cca:
//     R^-_(l), R^-_(l-1), ..., R_0, ..., R^+_(l-1), R^+_l
// cca-reverse:
//     R^+_(l), R^+_(l-1), ..., R_0, ..., R^-_(l-1), R^-_l
export type SphericalBasisOrder = 'gaussian' | 'cca' | 'cca-reverse';

export function normalizeBasicOrder(
    L: number,
    alpha: number[],
    order: SphericalBasisOrder
) {
    if (order === 'gaussian' || L === 0) return alpha;

    const ret: number[] = [alpha[L]];
    for (let l = 0; l < L; l++) {
        const a = alpha[L - l - 1],
            b = alpha[L + l + 1];
        if (order === 'cca') ret.push(b, a);
        else ret.push(a, b);
    }
    return ret;
}

export type SphericalFunc = (
    alpha: number[],
    x: number,
    y: number,
    z: number
) => number;

export const SphericalFunctions: SphericalFunc[] = [L0, L1, L2, L3, L4];

// L_i functions were auto-generated.

function L0(alpha: number[], x: number, y: number, z: number) {
    return alpha[0];
}

function L1(alpha: number[], x: number, y: number, z: number) {
    return alpha[0] * z + alpha[1] * x + alpha[2] * y;
}

function L2(alpha: number[], x: number, y: number, z: number) {
    const xx = x * x, yy = y * y, zz = z * z;
    return (
        alpha[0] * (-0.5 * xx - 0.5 * yy + zz) +
        alpha[1] * (1.7320508075688772 * x * z) +
        alpha[2] * (1.7320508075688772 * y * z) +
        alpha[3] * (0.8660254037844386 * xx - 0.8660254037844386 * yy) +
        alpha[4] * (1.7320508075688772 * x * y)
    );
}

function L3(alpha: number[], x: number, y: number, z: number) {
    const xx = x * x, yy = y * y, zz = z * z;
    const xxx = xx * x, yyy = yy * y, zzz = zz * z;
    return (
        alpha[0] * (-1.5 * xx * z - 1.5 * yy * z + zzz) +
        alpha[1] * (-0.6123724356957945 * xxx - 0.6123724356957945 * x * yy + 2.449489742783178 * x * zz) +
        alpha[2] * (-0.6123724356957945 * xx * y - 0.6123724356957945 * yyy + 2.449489742783178 * y * zz) +
        alpha[3] * (1.9364916731037085 * xx * z - 1.9364916731037085 * yy * z) +
        alpha[4] * (3.872983346207417 * x * y * z) +
        alpha[5] * (0.7905694150420949 * xxx - 2.3717082451262845 * x * yy) +
        alpha[6] * (2.3717082451262845 * xx * y - 0.7905694150420949 * yyy)
    );
}

function L4(alpha: number[], x: number, y: number, z: number) {
    const xx = x * x, yy = y * y, zz = z * z;
    const xxx = xx * x, yyy = yy * y, zzz = zz * z;
    const xxxx = xxx * x, yyyy = yyy * y, zzzz = zzz * z;
    return (
        alpha[0] * (0.375 * xxxx + 0.75 * xx * yy + 0.375 * yyyy - 3.0 * xx * zz - 3.0 * yy * zz + zzzz) +
        alpha[1] * (-2.3717082451262845 * xxx * z - 2.3717082451262845 * x * yy * z + 3.1622776601683795 * x * zzz) +
        alpha[2] * (-2.3717082451262845 * xx * y * z - 2.3717082451262845 * yyy * z + 3.1622776601683795 * y * zzz) +
        alpha[3] * (-0.5590169943749475 * xxxx + 0.5590169943749475 * yyyy + 3.3541019662496847 * xx * zz - 3.3541019662496847 * yy * zz) +
        alpha[4] * (-1.118033988749895 * xxx * y - 1.118033988749895 * x * yyy + 6.708203932499369 * x * y * zz) +
        alpha[5] * (2.091650066335189 * xxx * z + -6.274950199005566 * x * yy * z) +
        alpha[6] * (6.274950199005566 * xx * y * z + -2.091650066335189 * yyy * z) +
        alpha[7] * (0.739509972887452 * xxxx - 4.437059837324712 * xx * yy + 0.739509972887452 * yyyy) +
        alpha[8] * (2.958039891549808 * xxx * y + -2.958039891549808 * x * yyy)
    );
}
