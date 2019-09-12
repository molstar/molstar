/**
 * Copyright (c) 2018 Mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Camera } from '../camera';
import { Quat, Vec3 } from '../../mol-math/linear-algebra';
import { lerp } from '../../mol-math/interpolate';

export { CameraTransitionManager }

class CameraTransitionManager {
    private t = 0;

    private func: CameraTransitionManager.TransitionFunc = CameraTransitionManager.defaultTransition;
    private start = 0;
    inTransition = false;
    private durationMs = 0;
    private source: Camera.Snapshot = Camera.createDefaultSnapshot();
    private target: Camera.Snapshot = Camera.createDefaultSnapshot();
    private current = Camera.createDefaultSnapshot();

    apply(to: Partial<Camera.Snapshot>, durationMs: number = 0, transition?: CameraTransitionManager.TransitionFunc) {
        if (durationMs <= 0 || (typeof to.mode !== 'undefined' && to.mode !== this.camera.state.mode)) {
            this.finish(to);
            return;
        }

        Camera.copySnapshot(this.source, this.camera.state);
        Camera.copySnapshot(this.target, this.camera.state);
        Camera.copySnapshot(this.target, to);

        this.inTransition = true;
        this.func = transition || CameraTransitionManager.defaultTransition;
        this.start = this.t;
        this.durationMs = durationMs;
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
            this.finish(this.target!);
            return;
        }

        this.func(this.current, normalized, this.source, this.target);
        Camera.copySnapshot(this.camera.state, this.current);
    }

    constructor(private camera: Camera) {

    }
}

namespace CameraTransitionManager {
    export type TransitionFunc = (out: Camera.Snapshot, t: number, source: Camera.Snapshot, target: Camera.Snapshot) => void

    export function defaultTransition(out: Camera.Snapshot, t: number, source: Camera.Snapshot, target: Camera.Snapshot): void {
        Camera.copySnapshot(out, target);

        // Rotate direction
        const rot = Quat.identity();
        Quat.slerp(rot, rot, Quat.rotationTo(Quat.zero(), source.direction, target.direction), t);
        Vec3.transformQuat(out.direction, source.direction, rot);

        // Rotate up
        Quat.setIdentity(rot);
        Quat.slerp(rot, rot, Quat.rotationTo(Quat.zero(), source.up, target.up), t);
        Vec3.transformQuat(out.up, source.up, rot);

        // Lerp target
        Vec3.lerp(out.target, source.target, target.target, t);

        // Update position
        const dist = -lerp(Vec3.distance(source.position, source.target), Vec3.distance(target.position, target.target), t);
        Vec3.scale(out.position, out.direction, dist);
        Vec3.add(out.position, out.position, out.target);

        // Lerp other props
        out.zoom = lerp(source.zoom, target.zoom, t);
        out.fov = lerp(source.fov, target.fov, t);
        out.near = lerp(source.near, target.near, t);
        out.far = lerp(source.far, target.far, t);
        out.fogNear = lerp(source.fogNear, target.fogNear, t);
        out.fogFar = lerp(source.fogFar, target.fogFar, t);
    }
}