/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 *
 * adapted from https://github.com/d3/d3-ease
 */

const b1 = 4 / 11,
    b2 = 6 / 11,
    b3 = 8 / 11,
    b4 = 3 / 4,
    b5 = 9 / 11,
    b6 = 10 / 11,
    b7 = 15 / 16,
    b8 = 21 / 22,
    b9 = 63 / 64,
    b0 = 1 / b1 / b1;

export function bounceIn(t: number) {
    return 1 - bounceOut(1 - t);
}

export function bounceOut(t: number) {
    return (t = +t) < b1 ? b0 * t * t : t < b3 ? b0 * (t -= b2) * t + b4 : t < b6 ? b0 * (t -= b5) * t + b7 : b0 * (t -= b8) * t + b9;
}

export function bounceInOut(t: number) {
    return ((t *= 2) <= 1 ? 1 - bounceOut(1 - t) : bounceOut(t - 1) + 1) / 2;
}

//

export function circleIn(t: number) {
    return 1 - Math.sqrt(1 - t * t);
}

export function circleOut(t: number) {
    return Math.sqrt(1 - --t * t);
}

export function circleInOut(t: number) {
    return ((t *= 2) <= 1 ? 1 - Math.sqrt(1 - t * t) : Math.sqrt(1 - (t -= 2) * t) + 1) / 2;
}

//

export function cubicIn(t: number) {
    return t * t * t;
}

export function cubicOut(t: number) {
    return --t * t * t + 1;
}

export function cubicInOut(t: number) {
    return ((t *= 2) <= 1 ? t * t * t : (t -= 2) * t * t + 2) / 2;
}

//

export function expIn(t: number) {
    return Math.pow(2, 10 * t - 10);
}

export function expOut(t: number) {
    return 1 - Math.pow(2, -10 * t);
}

export function expInOut(t: number) {
    return ((t *= 2) <= 1 ? Math.pow(2, 10 * t - 10) : 2 - Math.pow(2, 10 - 10 * t)) / 2;
}

//

export function quadIn(t: number) {
    return t * t;
}

export function quadOut(t: number) {
    return t * (2 - t);
}

export function quadInOut(t: number) {
    return ((t *= 2) <= 1 ? t * t : --t * (2 - t) + 1) / 2;
}

//

const pi = Math.PI,
    halfPi = pi / 2;

export function sinIn(t: number) {
    return 1 - Math.cos(t * halfPi);
}

export function sinOut(t: number) {
    return Math.sin(t * halfPi);
}

export function sinInOut(t: number) {
    return (1 - Math.cos(pi * t)) / 2;
}

//
