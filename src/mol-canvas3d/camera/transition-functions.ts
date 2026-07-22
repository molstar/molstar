/**
 * Copyright (c) 2018-2026 Mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Adam Midlik <midlik@gmail.com>
 */

import { lerp } from '../../mol-math/interpolate';
import { Quat, Vec3 } from '../../mol-math/linear-algebra';
import { ParamDefinition as PD } from '../../mol-util/param-definition';
import { Camera } from '../camera';
import { CameraTransitionManager } from './transition';


/** Amount artificially added to the visible sphere radius in calculations for constant relative speed (dtarget / dt = (r + CONST_REL_SPEED_OFFSET) * const),
 * to avoid issues related to zero radius (1/0, log(0)) */
const CONST_REL_SPEED_OFFSET = 1;


/** Simple linear transition with constant absolute speed. */
export function transition_linear(out: Camera.Snapshot, t: number, source: Camera.Snapshot, target: Camera.Snapshot): void {
    return transition_linear_internal(out, t, t, source, target);
}

/** Linear transition with constant speed relative to visible sphere radius (move slower where zoomed-in more, dtarget / dt = (r + CONST_REL_SPEED_OFFSET) * const).
 * Rotational component of the transition has constant speed. */
export function transition_linear_relative(out: Camera.Snapshot, t: number, source: Camera.Snapshot, target: Camera.Snapshot): void {
    const rVisSource = visibleSphereRadius(source);
    const rVisTarget = visibleSphereRadius(target);
    const tTrans = constRelSpeedQuotientAdj_linRadIntp(rVisSource + CONST_REL_SPEED_OFFSET, rVisTarget + CONST_REL_SPEED_OFFSET, t);
    return transition_linear_internal(out, tTrans, t, source, target);
}

/** Linear transition allowing different transition quotient for distances (`tTrans`) and angles (`tRot`). */
function transition_linear_internal(out: Camera.Snapshot, tTrans: number, tRot: number, source: Camera.Snapshot, target: Camera.Snapshot): void {
    Camera.copySnapshot(out, target);

    // Interpolate fov & fog (use tRot as these are scale-independent)
    out.fov = lerp(source.fov, target.fov, tRot);
    out.fog = lerp(source.fog, target.fog, tRot);

    // Interpolate target
    Vec3.lerp(out.target, source.target, target.target, tTrans);

    // Interpolate distance (indirectly via visible sphere radius, for a more natural transition when FOV changes)
    const rVisSource = visibleSphereRadius(source);
    const rVisTarget = visibleSphereRadius(target);
    const rVis = lerp(rVisSource, rVisTarget, tTrans);
    const dist = cameraTargetDistance(rVis, out.mode, out.fov);

    // Interpolate direction and up
    interpolateCameraRotation(out, tRot, dist, source, target);

    // Interpolate radius
    out.radius = lerp(source.radius, target.radius, tTrans);
    out.radiusMax = lerp(source.radiusMax, target.radiusMax, tTrans);
}

const _sourceDirection = Vec3();
const _targetDirection = Vec3();
const _rotUp = Quat.identity();
const _rotDist = Quat.identity();

/** Interpolate camera direction and up, set camera distance from its target to `dist`. */
function interpolateCameraRotation(out: Camera.Snapshot, t: number, dist: number, source: Camera.Snapshot, target: Camera.Snapshot): void {
    // Rotate up
    Quat.rotationTo(_rotUp, source.up, target.up);
    Quat.slerp(_rotUp, Quat.Identity, _rotUp, t);
    Vec3.transformQuat(out.up, source.up, _rotUp);

    // Rotate between source and target direction
    Vec3.sub(_sourceDirection, source.position, source.target);
    Vec3.normalize(_sourceDirection, _sourceDirection);

    Vec3.sub(_targetDirection, target.position, target.target);
    Vec3.normalize(_targetDirection, _targetDirection);

    Quat.rotationTo(_rotDist, _sourceDirection, _targetDirection);
    Quat.slerp(_rotDist, Quat.Identity, _rotDist, t);
    Vec3.transformQuat(_sourceDirection, _sourceDirection, _rotDist);

    Vec3.scale(_sourceDirection, _sourceDirection, dist);
    Vec3.add(out.position, out.target, _sourceDirection);
}

/** "Leaping" camera transition with constant absolute speed.
 * Visible sphere radius increases during the transition so that for some t (0<=t<=1) the interpolated sphere will contain both source and target spheres.
 * When source and target spheres overlap, the transition is partially or fully linearized to avoid disturbing zoom-out.
 * Camera target position and rotational component of the transition use linear interpolation. */
export function transition_leap(out: Camera.Snapshot, t: number, source: Camera.Snapshot, target: Camera.Snapshot): void {
    Camera.copySnapshot(out, target);

    // Interpolate fov & fog
    out.fov = lerp(source.fov, target.fov, t);
    out.fog = lerp(source.fog, target.fog, t);

    // Interpolate distance (indirectly via visible sphere radius)
    const shift = Vec3.distance(source.target, target.target);
    const rVisSource = visibleSphereRadius(source);
    const rVisTarget = visibleSphereRadius(target);
    const rVis = leapingRadiusInterpolationSmart(rVisSource, rVisTarget, shift, t);
    const dist = cameraTargetDistance(rVis, out.mode, out.fov);

    // Interpolate target
    Vec3.lerp(out.target, source.target, target.target, t);

    // Interpolate direction and up
    interpolateCameraRotation(out, t, dist, source, target);

    // Interpolate radius
    out.radius = leapingRadiusInterpolationSmart(source.radius, target.radius, shift, t);
    out.radiusMax = leapingRadiusInterpolationSmart(source.radiusMax, target.radiusMax, shift, t);
}

/** "Leaping" camera transition with constant speed relative to visible sphere radius (move slower where zoomed-in more, dtarget / dt = (r + CONST_REL_SPEED_OFFSET) * const).
 * Rotational component of the transition uses linear interpolation. */
function transition_leap_relative(out: Camera.Snapshot, t: number, source: Camera.Snapshot, target: Camera.Snapshot): void {
    Camera.copySnapshot(out, target);

    // Interpolate fov & fog
    out.fov = lerp(source.fov, target.fov, t);
    out.fog = lerp(source.fog, target.fog, t);

    // Interpolate distance (indirectly via visible sphere radius)
    const shift = Vec3.distance(source.target, target.target);
    const rVisSource = visibleSphereRadius(source);
    const rVisTarget = visibleSphereRadius(target);
    const { r: rVis, q: tTrans } = getRadiusAndQuotientWithOffset(rVisSource, rVisTarget, shift, t);
    const dist = cameraTargetDistance(rVis, out.mode, out.fov);
    const rCorrection = rVis / lerp(rVisSource, rVisTarget, t);

    // Interpolate target
    Vec3.lerp(out.target, source.target, target.target, tTrans);

    // Interpolate direction and up
    interpolateCameraRotation(out, t, dist, source, target);

    // Interpolate radius
    out.radius = lerp(source.radius, target.radius, t) * rCorrection;
    out.radiusMax = lerp(source.radiusMax, target.radiusMax, t) * rCorrection;
}


export const TransitionFunctions = {
    'linear': transition_linear,
    'linear-relative': transition_linear_relative,
    'leap': transition_leap,
    'leap-relative': transition_leap_relative,
} satisfies Record<string, CameraTransitionManager.TransitionFunc>;

export type TransitionShape = keyof typeof TransitionFunctions;

export function getTransitionFn(shape: TransitionShape | undefined): CameraTransitionManager.TransitionFunc {
    if (!shape) return TransitionFunctions.linear;
    return TransitionFunctions[shape] ?? TransitionFunctions.linear;
}

export function TransitionShapeParamDefinition(defaultValue: TransitionShape): PD.Select<TransitionShape> {
    return PD.Select(
        defaultValue,
        Object.keys(TransitionFunctions).map(key => [key as TransitionShape, key]),
        { description: 'Camera transition trajectory shape. "linear": interpolates along a straight line with constant absolute speed; "linear-relative": like "linear" but moves relatively slower when zoomed-in more; "leap": zooms out during the transition to capture both the initial and the final camera target (becomes linear if the targets are near); "leap-relative": like "leap" but moves relatively slower when zoomed-in more.' }
    );
}


/** Sphere radius "interpolation" method which increases the radius during transition so that for some t (0<=t<=1) the interpolated sphere will contain both source and target spheres.
 * Assumes that target itself uses linear interpolation.
 * `rA`, `rB` - radius of source (t=0) and target (t=1) sphere;
 * `dist` - distance between centers of source and target sphere. */
function leapingRadiusInterpolationCubic(rA: number, rB: number, dist: number, t: number): number {
    if (dist === 0) {
        return lerp(rA, rB, t);
    }
    if (rB >= dist + rA) { // Sphere 1 fully contains sphere 0
        const alpha = dist / (rB - rA);
        return lerp(rA, rB, niceCubic(t, alpha));
    }
    if (rA >= dist + rB) { // Sphere 0 fully contains sphere 1
        const alpha = dist / (rA - rB);
        return lerp(rB, rA, niceCubic(1 - t, alpha));
    }
    const tmax = (dist - rA + rB) / 2 / dist;
    const rmax = (dist + rA + rB) / 2;
    if (t <= tmax) {
        return lerp(rA, rmax, niceCubic(t / tmax));
    } else {
        return lerp(rB, rmax, niceCubic((1 - t) / (1 - tmax)));
    }
}

/** Sphere radius "interpolation" method similar to `leapingRadiusInterpolationCubic`,
 * but the radius increases less when the source and target spheres overlap, and becomes linear when at least one of the spheres contains the center of the other
 * (this is to avoid disturbing zoom-out when the source and target are near). */
function leapingRadiusInterpolationSmart(rA: number, rB: number, dist: number, t: number): number {
    const overlapFactor = relativeSphereOverlap(rA, rB, dist);
    if (overlapFactor <= 0) return leapingRadiusInterpolationCubic(rA, rB, dist, t); // spheres not overlapping
    if (overlapFactor >= 1) return lerp(rA, rB, t); // one of the spheres contains the center of the other
    return lerp(leapingRadiusInterpolationCubic(rA, rB, dist, t), lerp(rA, rB, t), overlapFactor);
}

/** Arbitrary measure of how much two spheres overlap (>0 when spheres do not overlap, >=1 when at least of the spheres contains the center of the other) */
function relativeSphereOverlap(rA: number, rB: number, dist: number): number {
    const overlap = rA + rB - dist;
    if (rA === 0 || rB === 0) {
        return overlap >= 0 ? Infinity : -Infinity;
    }
    return overlap / Math.min(rA, rB);
}

/** Auxiliary cubic function that goes from y(0)=0 to y(1)=1.
 * When alpha=1, it is a curve with inflection point in 0 and stationary point in 1.
 * When alpha=0, it becomes a linear function. */
function niceCubic(x: number, alpha: number = 1) {
    return (1 + 0.5 * alpha) * x - 0.5 * alpha * x ** 3;
}

/** Return the radius of the largest sphere centered in camera.target which is fully in FOV */
function visibleSphereRadius(camera: Camera.Snapshot): number {
    const distance = Vec3.distance(camera.target, camera.position);
    if (camera.mode === 'orthographic')
        return distance * Math.tan(camera.fov / 2);
    else // perspective
        return distance * Math.sin(camera.fov / 2);
}
/** Return the distance of camera from the center of a sphere with radius `visRadius` so that the sphere just fits into FOV */
function cameraTargetDistance(visRadius: number, mode: Camera.Mode, fov: number): number {
    if (mode === 'orthographic')
        return visRadius / Math.tan(fov / 2);
    else // perspective
        return visRadius / Math.sin(fov / 2);
}

/** This adjustment to transition quotient (0-1) ensures that transition speed relative to radius is constant during the transition.
 * Only to be applied to position and radius interpolation, not needed for angles.
 * This function assumes linear interpolation of radius. For other interpolation methods, more complicated formula is needed. */
function constRelSpeedQuotientAdj_linRadIntp(rA: number, rB: number, t: number): number {
    if (isZero((rA - rB) / (rA + rB))) {
        // Special case for rA===rB
        return t;
    }
    return rA / (rB - rA) * ((rB / rA) ** t - 1); // = a / b * ((1 + b / a) ** t - 1), where a = rA, b = rB - rA
}

function isZero(x: number) {
    return Math.abs(x) <= 0.001;
}

/** Exponential function `r * Math.exp(k * x)` with added "hump" so that derivation in x=1 is increased (alpha>0) or decreased (alpha<0).  */
function fx(r: number, k: number, alpha: number, x: number) {
    return r * Math.exp(k * x) * (1 + alpha / 2 * x ** 3 - alpha / 2 * x);
}

/** Integral of `fx` */
function integralFx(r: number, k: number, alpha: number, x: number) {
    return r * Math.exp(k * x) * (
        alpha / (2 * k) * x ** 3
        - 3 * alpha / (2 * k ** 2) * x ** 2
        + (3 * alpha / k ** 3 - alpha / (2 * k)) * x
        + 1 / k + alpha / (2 * k ** 2) - 3 * alpha / k ** 4);
}

/** Compute visible sphere radius `r` and transition quotient `q` (for camera target interpolation) at time `t` (0<=t<=1)
 * for leap-relative transition. Avoid issues with zero radius by adding `CONST_REL_SPEED_OFFSET`. */
function getRadiusAndQuotientWithOffset(rA: number, rB: number, dist: number, t: number) {
    const { r, q } = getRadiusAndQuotient(rA + CONST_REL_SPEED_OFFSET, rB + CONST_REL_SPEED_OFFSET, dist, t);
    return { r: r - CONST_REL_SPEED_OFFSET, q };
}

/** Compute visible sphere radius `r` and transition quotient `q` (for camera target interpolation) at time `t` (0<=t<=1)
 * for leap-relative transition, without adding radius offset. Assumes initial and final radii (`rA`, `rB`) are non-zero. */
function getRadiusAndQuotient(rA: number, rB: number, dist: number, t: number) {
    const overlap = relativeSphereOverlap(rA, rB, dist);
    /** Non-linearity coefficient (0 -> purely linear transition, 1 -> purely leaping transition) */
    const beta = Math.max(0, Math.min(1, 1 - overlap));
    if (beta === 0) {
        // Linear
        if (isZero(rA - rB)) {
            return { r: lerp(rA, rB, t), q: t };
        } else {
            const log_rA = Math.log(rA);
            const log_rB = Math.log(rB);
            const k0 = log_rB - log_rA;
            const alpha = 0;
            const r = fx(rA, k0, alpha, t);
            const F0 = integralFx(rA, k0, alpha, 0);
            const F1 = integralFx(rA, k0, alpha, 1);
            const c = 1 / (F1 - F0);
            const q = c * (integralFx(rA, k0, alpha, t) - F0);
            return { r, q };
        }
    } else {
        // Leaping
        /** Maximum radius if the transition were fully leaping (beta=1) */
        const rC_ideal = (dist + rA + rB) / 2;
        const log_rA = Math.log(rA);
        const log_rB = Math.log(rB);
        const log_rC_ideal = Math.log(rC_ideal);
        /** Time where maximum radius happens */
        const tau = (log_rC_ideal - log_rA) / ((log_rC_ideal - log_rA) + (log_rC_ideal - log_rB));
        /** Real maximum radius (partially linearized) */
        const log_rC = lerp(lerp(log_rA, log_rB, tau), log_rC_ideal, beta);
        const k0 = log_rB - log_rA;
        const kA = log_rC - log_rA;
        const kB = log_rC - log_rB;
        const alphaA = k0 * (1 - beta) * tau - kA;
        const alphaB = -k0 * (1 - beta) * (1 - tau) - kB;

        const FA0 = integralFx(rA, kA, alphaA, 0);
        const FA1 = integralFx(rA, kA, alphaA, 1);
        const FB0 = integralFx(rB, kB, alphaB, 0);
        const FB1 = integralFx(rB, kB, alphaB, 1);
        const c = 1 / (tau * (FA1 - FA0) + (1 - tau) * (FB1 - FB0));
        if (t <= tau) {
            const x = t / tau;
            const r = fx(rA, kA, alphaA, x);
            const q = c * tau * (integralFx(rA, kA, alphaA, x) - FA0);
            return { r, q };
        } else {
            const x = (1 - t) / (1 - tau);
            const r = fx(rB, kB, alphaB, x);
            const q = 1 - c * (1 - tau) * (integralFx(rB, kB, alphaB, x) - FB0);
            return { r, q };
        }
    }
}
