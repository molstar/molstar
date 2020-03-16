/**
 * Copyright (c) 2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Sphere3D } from '../../mol-math/geometry';
import { StructureElement } from '../../mol-model/structure';
import { PluginContext } from '../../mol-plugin/context';
import { PrincipalAxes } from '../../mol-math/linear-algebra/matrix/principal-axes';
import { Camera } from '../../mol-canvas3d/camera';

// TODO: make this customizable somewhere?
const DefaultCameraFocusOptions = {
    minRadius: 5,
    extraRadius: 6,
    durationMs: 250
}

export type CameraFocusOptions = typeof DefaultCameraFocusOptions

export class CameraManager {
    focusLoci(loci: StructureElement.Loci, options?: Partial<CameraFocusOptions>) {
        // TODO: allow computation of principal axes here?
        // perhaps have an optimized function, that does exact axes small Loci and approximate/sampled from big ones?

        const { extraRadius, minRadius, durationMs } = { ...DefaultCameraFocusOptions, ...options };
        const { sphere } = StructureElement.Loci.getBoundary(loci);
        const radius = Math.max(sphere.radius + extraRadius, minRadius);
        this.plugin.canvas3d?.camera.focus(sphere.center, radius, durationMs);
    }

    focusSphere(sphere: Sphere3D, options?: Partial<CameraFocusOptions> & { principalAxes?: PrincipalAxes }) {
        const { extraRadius, minRadius, durationMs } = { ...DefaultCameraFocusOptions, ...options };
        const radius = Math.max(sphere.radius + extraRadius, minRadius);

        if (options?.principalAxes) {
            const { origin, dirA, dirC } = options?.principalAxes.boxAxes;
            this.plugin.canvas3d?.camera.focus(origin, radius, durationMs, dirA, dirC);
        } else {
            this.plugin.canvas3d?.camera.focus(sphere.center, radius, durationMs);
        }
    }

    setSnapshot(snapshot: Partial<Camera.Snapshot>, durationMs?: number) {
        // TODO: setState and requestCameraReset are very similar now: unify them?
        this.plugin.canvas3d?.camera.setState(snapshot, durationMs);
    }

    reset(snapshot?: Partial<Camera.Snapshot>, durationMs?: number) {
        this.plugin.canvas3d?.requestCameraReset({ snapshot, durationMs });
    }

    constructor(readonly plugin: PluginContext) {
    }
}