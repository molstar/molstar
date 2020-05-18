/**
 * Copyright (c) 2019-2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Sphere3D } from '../../mol-math/geometry';
import { PluginContext } from '../../mol-plugin/context';
import { PrincipalAxes } from '../../mol-math/linear-algebra/matrix/principal-axes';
import { Camera } from '../../mol-canvas3d/camera';
import { Loci } from '../../mol-model/loci';
import { BoundaryHelper } from '../../mol-math/geometry/boundary-helper';
import { GraphicsRenderObject } from '../../mol-gl/render-object';
import { StructureElement } from '../../mol-model/structure';

// TODO: make this customizable somewhere?
const DefaultCameraFocusOptions = {
    minRadius: 5,
    extraRadius: 6,
    durationMs: 250
};

export type CameraFocusOptions = typeof DefaultCameraFocusOptions

export class CameraManager {
    private boundaryHelper = new BoundaryHelper('98');

    private transformedLoci(loci: Loci) {
        if (StructureElement.Loci.is(loci)) {
            // use decorated (including 3d transforms) parent structure
            const parent = this.plugin.helpers.substructureParent.get(loci.structure)?.obj?.data;
            if (parent) loci = StructureElement.Loci.remap(loci, parent);
        }
        return loci;
    }

    focusRenderObjects(objects?: ReadonlyArray<GraphicsRenderObject>, options?: Partial<CameraFocusOptions>) {
        if (!objects) return;
        const spheres: Sphere3D[] = [];

        for (const o of objects) {
            const s = o.values.boundingSphere.ref.value;
            if (s.radius === 0) continue;
            spheres.push(s);
        }

        this.focusSpheres(spheres, s => s, options);
    }

    focusLoci(loci: Loci | Loci[], options?: Partial<CameraFocusOptions>) {
        // TODO: allow computation of principal axes here?
        // perhaps have an optimized function, that does exact axes small Loci and approximate/sampled from big ones?

        let sphere: Sphere3D | undefined;

        if (Array.isArray(loci) && loci.length > 1) {
            const spheres = [];
            for (const l of loci) {
                const s = Loci.getBoundingSphere(this.transformedLoci(l));
                if (s) spheres.push(s);
            }

            if (spheres.length === 0) return;

            this.boundaryHelper.reset();
            for (const s of spheres) {
                this.boundaryHelper.includeSphere(s);
            }
            this.boundaryHelper.finishedIncludeStep();
            for (const s of spheres) {
                this.boundaryHelper.radiusSphere(s);
            }
            sphere = this.boundaryHelper.getSphere();
        } else if (Array.isArray(loci)) {
            if (loci.length === 0) return;
            sphere = Loci.getBoundingSphere(this.transformedLoci(loci[0]));
        } else {
            sphere = Loci.getBoundingSphere(this.transformedLoci(loci));
        }

        if (sphere) {
            this.focusSphere(sphere, options);
        }
    }

    focusSpheres<T>(xs: ReadonlyArray<T>, sphere: (t: T) => Sphere3D | undefined, options?: Partial<CameraFocusOptions>) {
        const spheres = [];

        for (const x of xs) {
            const s = sphere(x);
            if (s) spheres.push(s);
        }

        if (spheres.length === 0) return;
        if (spheres.length === 1) return this.focusSphere(spheres[0], options);

        this.boundaryHelper.reset();
        for (const s of spheres) {
            this.boundaryHelper.includeSphere(s);
        }
        this.boundaryHelper.finishedIncludeStep();
        for (const s of spheres) {
            this.boundaryHelper.radiusSphere(s);
        }
        this.focusSphere(this.boundaryHelper.getSphere(), options);
    }

    focusSphere(sphere: Sphere3D, options?: Partial<CameraFocusOptions> & { principalAxes?: PrincipalAxes }) {
        const { extraRadius, minRadius, durationMs } = { ...DefaultCameraFocusOptions, ...options };
        const radius = Math.max(sphere.radius + extraRadius, minRadius);

        if (options?.principalAxes) {
            const { origin, dirA, dirC } = options?.principalAxes.boxAxes;
            const snapshot = this.plugin.canvas3d?.camera.getFocus(origin, radius, dirA, dirC);
            this.plugin.canvas3d?.requestCameraReset({ durationMs, snapshot });
            // this.plugin.canvas3d?.camera.focus(origin, radius, durationMs, dirA, dirC);
        } else {
            const snapshot = this.plugin.canvas3d?.camera.getFocus(sphere.center, radius);
            this.plugin.canvas3d?.requestCameraReset({ durationMs, snapshot });

            // this.plugin.canvas3d?.camera.focus(sphere.center, radius, durationMs);
        }
    }

    setSnapshot(snapshot: Partial<Camera.Snapshot>, durationMs?: number) {
        // TODO: setState and requestCameraReset are very similar now: unify them?
        this.plugin.canvas3d?.requestCameraReset({ snapshot, durationMs });
    }

    reset(snapshot?: Partial<Camera.Snapshot>, durationMs?: number) {
        this.plugin.canvas3d?.requestCameraReset({ snapshot, durationMs });
    }

    constructor(readonly plugin: PluginContext) {
    }
}