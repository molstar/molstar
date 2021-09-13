/**
 * Copyright (c) 2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 *
 * mostly adapted from https://gist.github.com/imbcmdth/6338194
 * which is ported from https://code.google.com/archive/p/fastapprox/ (BSD licensed)
 */

const _a_fastPow2 = new ArrayBuffer(4);
const _i_fastPow2 = new Int32Array(_a_fastPow2);
const _f_fastPow2 = new Float32Array(_a_fastPow2);

export function fastPow2(v: number) {
    const offset = (v < 0) ? 1 : 0;
    const clipNumber = (v < -126) ? -126 : v;
    const w = clipNumber | 0;
    const z = clipNumber - w + offset;
    _i_fastPow2[0] = ((1 << 23) * (clipNumber + 121.2740575 + 27.7280233 / (4.84252568 - z) - 1.49012907 * z));
    return _f_fastPow2[0];
}

const _a_fasterPow2 = new ArrayBuffer(4);
const _i_fasterPow2 = new Int32Array(_a_fasterPow2);
const _f_fasterPow2 = new Float32Array(_a_fasterPow2);

export function fasterPow2(v: number) {
    const clipNumber = (v < -126) ? -126 : v;
    _i_fasterPow2[0] = ((1 << 23) * (clipNumber + 126.94269504));
    return _f_fasterPow2[0];
}

export function fastExp(v: number) {
    return fastPow2(1.442695040 * v);
}

export function fasterExp(v: number) {
    return fasterPow2(1.442695040 * v);
}

const _a_fastLog2 = new ArrayBuffer(8);
const _i_fastLog2 = new Int32Array(_a_fastLog2);
const _f_fastLog2 = new Float32Array(_a_fastLog2);

export function fastLog2(v: number) {
    _f_fastLog2[0] = v;
    _i_fastLog2[1] = (_i_fastLog2[0] & 0x007FFFFF) | 0x3f000000;
    const t = _i_fastLog2[0] * 1.1920928955078125e-7;
    return t - 124.22551499 - 1.498030302 * _f_fastLog2[1] - 1.72587999 / (0.3520887068 + _f_fastLog2[1]);
};

const _a_fasterLog2 = new ArrayBuffer(4);
const _i_fasterLog2 = new Int32Array(_a_fasterLog2);
const _f_fasterLog2 = new Float32Array(_a_fasterLog2);

export function fasterLog2(v: number) {
    _f_fasterLog2[0] = v;
    const t = _i_fasterLog2[0] * 1.1920928955078125e-7;
    return t - 126.94269504;
}

export function fastLog(v: number) {
    return 0.6931471805599453 * fastLog2(v);
}

export function fasterLog(v: number) {
    return 0.6931471805599453 * fasterLog2(v);
}

export function fastLog10(v: number) {
    return 0.30102999566398114 * fastLog2(v);
}

export function fasterLog10(v: number) {
    return 0.30102999566398114 * fasterLog2(v);
}

export function fastSinh(v: number) {
    return 0.5 * (fastExp(v) - fastExp(-v));
}

export function fasterSinh(v: number) {
    return 0.5 * (fasterExp(v) - fasterExp(-v));
}

export function fastCosh(v: number) {
    return 0.5 * (fastExp(v) + fastExp(-v));
}

export function fasterCosh(v: number) {
    return 0.5 * (fasterExp(v) + fasterExp(-v));
}

export function fastTanh(v: number) {
    return -1.0 + 2.0 / (1.0 + fastExp(-2.0 * v));
}

export function fasterTanh(v: number) {
    return -1.0 + 2.0 / (1.0 + fasterExp(-2.0 * v));
}

const halfPi = Math.PI / 2;
const twoPi = 2 * Math.PI;
const invTwoPi = 1 / (2 * Math.PI);
const twoOverPi = 2 / Math.PI;
const fourOverPi = 4 / Math.PI;
const fourOverPiSq = 4 / (Math.PI * Math.PI);
const halfPiMinusTwoPi = Math.PI / 2 - 2 * Math.PI;

const _q_fastHalfSin = 0.78444488374548933;
const _a_fastHalfSin = new ArrayBuffer(16);
const _i_fastHalfSin = new Int32Array(_a_fastHalfSin);
const _f_fastHalfSin = new Float32Array(_a_fastHalfSin);

function fastHalfSin(v: number) {
    _f_fastHalfSin[0] = 0.20363937680730309;
    _f_fastHalfSin[1] = 0.015124940802184233;
    _f_fastHalfSin[2] = -0.0032225901625579573;
    _f_fastHalfSin[3] = v;
    const sign = _i_fastHalfSin[3] & 0x80000000;
    _i_fastHalfSin[3] = _i_fastHalfSin[3] & 0x7FFFFFFF;
    const qpprox = fourOverPi * v - fourOverPiSq * v * _f_fastHalfSin[3];
    const qpproxsq = qpprox * qpprox;
    _i_fastHalfSin[0] |= sign;
    _i_fastHalfSin[1] |= sign;
    _i_fastHalfSin[2] ^= sign;
    return _q_fastHalfSin * qpprox + qpproxsq * (_f_fastHalfSin[0] + qpproxsq * (_f_fastHalfSin[1] + qpproxsq * _f_fastHalfSin[2]));
}

const _q_fasterHalfSin = 0.78444488374548933;
const _a_fasterHalfSin = new ArrayBuffer(8);
const _i_fasterHalfSin = new Int32Array(_a_fasterHalfSin);
const _f_fasterHalfSin = new Float32Array(_a_fasterHalfSin);

function fasterHalfSin(v: number) {
    _f_fasterHalfSin[0] = 0.22308510060189463;
    _f_fasterHalfSin[1] = v;
    const sign = _i_fasterHalfSin[1] & 0x80000000;
    _i_fasterHalfSin[1] &= 0x7FFFFFFF;
    const qpprox = fourOverPi * v - fourOverPiSq * v * _f_fasterHalfSin[1];
    _i_fasterHalfSin[0] |= sign;
    return qpprox * (_q_fasterHalfSin + _f_fasterHalfSin[0] * qpprox);
}

export function fastSin(v: number) {
    const k = (v * invTwoPi) | 0;
    const half = (v < 0) ? -0.5 : 0.5;
    return fastHalfSin((half + k) * twoPi - v);
}

export function fasterSin(v: number) {
    const k = (v * invTwoPi) | 0;
    const half = (v < 0) ? -0.5 : 0.5;
    return fasterHalfSin((half + k) * twoPi - v);
}

export function fastCos(v: number) {
    return fastSin(v + halfPi);
}

export function fasterCos(v: number) {
    return fasterSin(v + halfPi);
}

function fastHalfCos(v: number) {
    const offset = (v > halfPi) ? halfPiMinusTwoPi : halfPi;
    return fastHalfSin(v + offset);
}

const _p_fasterHalfCos = 0.54641335845679634;
const _a_fasterHalfCos = new ArrayBuffer(4);
const _i_fasterHalfCos = new Int32Array(_a_fasterHalfCos);
const _f_fasterHalfCos = new Float32Array(_a_fasterHalfCos);

function fasterHalfCos(v: number) {
    _f_fasterHalfCos[0] = v;
    _i_fasterHalfCos[0] &= 0x7FFFFFFF;
    const qpprox = 1.0 - twoOverPi * _f_fasterHalfCos[0];
    return qpprox + _p_fasterHalfCos * qpprox * (1.0 - qpprox * qpprox);
}

export function fastTan(v: number) {
    const k = (v * invTwoPi) | 0;
    const half = (v < 0) ? -0.5 : 0.5;
    const x = v - (half + k) * twoPi;
    return fastHalfSin(x) / fastHalfCos(x);
}

export function fasterTan(v: number) {
    const k = (v * invTwoPi) | 0;
    const half = (v < 0) ? -0.5 : 0.5;
    const x = v - (half + k) * twoPi;
    return fasterHalfSin(x) / fasterHalfCos(x);
}

const piOverFour = Math.PI / 4;

/**
 * Adapted from:
 * "Efficient approximations for the arctangent function"
 * Rajan, S. Sichun Wang Inkol, R. Joyal, A., May 2006
 */
export function fastAtan(v: number) {
    return piOverFour * v - v * (Math.abs(v) - 1) * (0.2447 + 0.0663 * Math.abs(v));
}

export function fastAtan2(y: number, x: number) {
    // reduce range to [-1, 1] by flipping y/x so the larger is up
    let t = Math.abs(x); // used to undo flipping
    let opposite = Math.abs(y);
    const adjacent = Math.max(t, opposite);
    opposite = Math.min(t, opposite);

    t = fastAtan(opposite / adjacent);
    // undo flipping
    t = Math.abs(y) > Math.abs(x) ? halfPi - t : t;
    t = x < 0.0 ? Math.PI - t : t;
    t = y < 0.0 ? -t : t;
    return t;
}