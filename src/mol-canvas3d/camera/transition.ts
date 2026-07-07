/**
 * Copyright (c) 2018-2026 Mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Adam Midlik <midlik@gmail.com>
 */

import { EasingFunction, getEasingFn } from '../../mol-math/easing';
import { Camera } from '../camera';
import { getTransitionFn, TransitionShape } from './transition-functions';

export { CameraTransitionManager };

export interface CameraTransitionOptions {
    /** If present, approximates the transion between, [current] -> [keyframes] -> -> [target] */
    keyframes?: CameraTransitionManager.TransitionKeyframes,
    /** Global easing, if easing is specified for keyframes, the "end" frame value is used  */
    easing?: EasingFunction,
    /** Defines shape of the camera trajectory during transition */
    shape?: TransitionShape,
}

class CameraTransitionManager {
    private t = 0;

    private func: CameraTransitionManager.TransitionFunc = CameraTransitionManager.defaultTransition;
    private start = 0;
    inTransition = false;
    private durationMs = 0;
    private _source: Camera.Snapshot = Camera.createDefaultSnapshot();
    private _target: Camera.Snapshot = Camera.createDefaultSnapshot();
    private _options: CameraTransitionOptions | undefined = void 0;
    private _current = Camera.createDefaultSnapshot();

    get source(): Readonly<Camera.Snapshot> { return this._source; }
    get target(): Readonly<Camera.Snapshot> { return this._target; }

    apply(
        to: Partial<Camera.Snapshot>,
        durationMs: number = 0,
        transition?: CameraTransitionManager.TransitionFunc,
        options?: CameraTransitionOptions,
    ) {
        if (!this.inTransition || durationMs > 0) {
            Camera.copySnapshot(this._source, this.camera.state);
        }

        if (!this.inTransition) {
            Camera.copySnapshot(this._target, this.camera.state);
        }

        Camera.copySnapshot(this._target, to);

        if (this._target.radius > this._target.radiusMax) {
            this._target.radius = this._target.radiusMax;
        }

        if (this._target.radius < 0.01) this._target.radius = 0.01;
        if (this._target.radiusMax < 0.01) this._target.radiusMax = 0.01;

        if (!this.inTransition && durationMs <= 0 || (typeof to.mode !== 'undefined' && to.mode !== this.camera.state.mode)) {
            this.finish(this._target);
            return;
        }

        this.inTransition = true;
        this.func = transition || CameraTransitionManager.defaultTransition;
        this._options = options;

        if (!this.inTransition || durationMs > 0) {
            this.start = this.t;
            this.durationMs = durationMs;
        }
    }

    tick(t: number) {
        this.t = t;
        this.update();
    }

    private finish(to: Partial<Camera.Snapshot>) {
        Camera.copySnapshot(this.camera.state, to);
        this.inTransition = false;
    }

    private update() {
        if (!this.inTransition) return;

        const normalized = Math.min((this.t - this.start) / this.durationMs, 1);
        if (normalized === 1) {
            this.finish(this._target!);
            return;
        }

        this.func(this._current, normalized, this._source, this._target, this._options);
        Camera.copySnapshot(this.camera.state, this._current);
    }

    constructor(private camera: Camera) {

    }
}

namespace CameraTransitionManager {
    export type TransitionKeyframes = { t: number, snapshot: Partial<Camera.Snapshot>, easing?: EasingFunction, shape?: TransitionShape }[]
    export type TransitionFunc = (out: Camera.Snapshot, t: number, source: Camera.Snapshot, target: Camera.Snapshot, options?: CameraTransitionOptions) => void

    let _tempSource: Camera.Snapshot | undefined = undefined;
    let _tempTarget: Camera.Snapshot | undefined = undefined;

    export function defaultTransition(
        out: Camera.Snapshot,
        t_: number,
        source_: Camera.Snapshot,
        target_: Camera.Snapshot,
        options?: CameraTransitionOptions
    ): void {
        let sourcePartial: Partial<Camera.Snapshot> = source_;
        let targetPartial: Partial<Camera.Snapshot> = target_;

        let tStart = 0;
        let tEnd = 1;
        let easingKind = options?.easing;
        let shapeKind = options?.shape;

        const keyframes = options?.keyframes;
        if (keyframes && keyframes.length > 0) {
            for (let i = 0; i < keyframes.length; i++) {
                const keyframe = keyframes[i];
                if (t_ >= keyframe.t) {
                    sourcePartial = keyframe.snapshot;
                    tStart = keyframe.t;
                    break;
                }
            }
            for (let i = 0; i < keyframes.length; i++) {
                const keyframe = keyframes[i];
                if (keyframe.t >= t_) {
                    targetPartial = keyframe.snapshot;
                    tEnd = keyframe.t;
                    easingKind = keyframe.easing ?? easingKind;
                    shapeKind = keyframe.shape ?? shapeKind;
                    break;
                }
            }
        }

        const easing = getEasingFn(easingKind);
        const t = easing((t_ - tStart) / (tEnd - tStart));

        if (!_tempSource) _tempSource = Camera.createDefaultSnapshot();
        if (!_tempTarget) _tempTarget = Camera.createDefaultSnapshot();

        Camera.copySnapshot(_tempSource, source_);
        Camera.copySnapshot(_tempSource, sourcePartial);
        Camera.copySnapshot(_tempTarget, target_);
        Camera.copySnapshot(_tempTarget, targetPartial);

        const transition = getTransitionFn(shapeKind);
        transition(out, t, _tempSource, _tempTarget);
    }
}