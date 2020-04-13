/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

export function normalize (value: number, min: number, max: number) {
    return (value - min) / (max - min);
}

export function clamp (value: number, min: number, max: number) {
    return Math.max(min, Math.min(max, value));
}

export function pclamp (value: number) {
    return clamp(value, 0, 100);
}

export function saturate (value: number) {
    return clamp(value, 0, 1);
}

export function damp (value: number, dampingFactor: number) {
    const dampedValue = value * dampingFactor;
    return Math.abs(dampedValue) < 0.1 ? 0 : dampedValue;
}

export function lerp (start: number, stop: number, alpha: number) {
    return start + (stop - start) * alpha;
}

/** Catmul-Rom spline */
export function spline (p0: number, p1: number, p2: number, p3: number, t: number, tension: number) {
    const v0 = (p2 - p0) * tension;
    const v1 = (p3 - p1) * tension;
    const t2 = t * t;
    const t3 = t * t2;
    return (2 * p1 - 2 * p2 + v0 + v1) * t3 + (-3 * p1 + 3 * p2 - 2 * v0 - v1) * t2 + v0 * t + p1;
}

export function quadraticBezier(p0: number, p1: number, p2: number, t: number) {
    const k = 1 - t;
    return (k * k * p0) + (2 * k * t * p1) + (t * t * p2);
}

export function smoothstep (min: number, max: number, x: number) {
    x = saturate(normalize(x, min, max));
    return x * x * (3 - 2 * x);
}

export function smootherstep (min: number, max: number, x: number) {
    x = saturate(normalize(x, min, max));
    return x * x * x * (x * (x * 6 - 15) + 10);
}

export function smootheststep (min: number, max: number, x: number) {
    x = saturate(normalize(x, min, max));
    return -20 * Math.pow(x, 7) + 70 * Math.pow(x, 6) - 84 * Math.pow(x, 5) + 35 * Math.pow(x, 4);
}

export function almostIdentity (value: number, start: number, stop: number) {
    if (value > start) return value;
    const a = 2 * stop - start;
    const b = 2 * start - 3 * stop;
    const t = value / start;
    return (a * t + b) * t * t + stop;
}