/**
 * Copyright (c) 2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 *
 * adapted from https://gist.github.com/imbcmdth/6338194
 * which is ported from https://code.google.com/archive/p/fastapprox/ (BSD licensed)
 */

const _a_fasterPow2 = new ArrayBuffer(4);
const _i_fasterPow2 = new Int32Array(_a_fasterPow2);
const _f_fasterPow2 = new Float32Array(_a_fasterPow2);

export function fasterPow2(v: number) {
    const clipNumber = (v < -126) ? -126 : v;
    _i_fasterPow2[0] = ((1 << 23) * (clipNumber + 126.94269504));
    return _f_fasterPow2[0];
}

export function fasterExp(v: number) {
    return fasterPow2(1.442695040 * v);
}

const _a_fasterLog2 = new ArrayBuffer(4);
const _i_fasterLog2 = new Int32Array(_a_fasterLog2);
const _f_fasterLog2 = new Float32Array(_a_fasterLog2);

export function fasterLog2(v: number) {
    _f_fasterLog2[0] = v;
    const t = _i_fasterLog2[0] * 1.1920928955078125e-7;
    return t - 126.94269504;
}

export function fasterLog(v: number) {
    return 0.6931471805599453 * fasterLog2(v);
}

export function fasterLog10(v: number) {
    return 0.30102999566398114 * fasterLog2(v);
}