/**
 * Copyright (c) 2019-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author Ke Ma <mark.ma@rcsb.org>
 * @author Adam Midlik <midlik@gmail.com>
 */

import { Camera } from '../../mol-canvas3d/camera';
import { GraphicsRenderObject } from '../../mol-gl/render-object';
import { Sphere3D } from '../../mol-math/geometry';
import { BoundaryHelper } from '../../mol-math/geometry/boundary-helper';
import { Mat3 } from '../../mol-math/linear-algebra';
import { leastObstructedDirection } from '../../mol-math/linear-algebra/3d/optimize-direction';
import { Vec3 } from '../../mol-math/linear-algebra/3d/vec3';
import { PrincipalAxes } from '../../mol-math/linear-algebra/matrix/principal-axes';
import { Loci } from '../../mol-model/loci';
import { Structure, StructureElement, StructureProperties } from '../../mol-model/structure';
import { PluginContext } from '../../mol-plugin/context';
import { PluginState } from '../../mol-plugin/state';
import { PluginStateObject } from '../objects';
import { pcaFocus } from './focus-camera/focus-first-residue';
import { getFocusSnapshot } from './focus-camera/focus-object';
import { changeCameraRotation, structureLayingTransform } from './focus-camera/orient-axes';

// TODO: make this customizable somewhere?
const DefaultCameraFocusOptions = {
    minRadius: 1,
    extraRadius: 4,
    durationMs: 250,
};

export type CameraFocusOptions = typeof DefaultCameraFocusOptions
export type CameraFocusLociOptions = CameraFocusOptions & { optimizeDirection?: boolean, optimizeDirectionUp?: Vec3 }

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

    private getFocusSphere(loci: Loci | Loci[]) {
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

        return sphere;
    }

    private focusLociOptimized(loci: Loci | Loci[], options?: Partial<CameraFocusOptions & { up?: Vec3 }>) {
        const { canvas3d } = this.plugin;
        if (!canvas3d) return;

        const sphere = this.getFocusSphere(loci);
        if (!sphere) return;

        const lociArray = Array.isArray(loci) ? loci : [loci];
        const positions: { x: number[], y: number[], z: number[] } = { x: [], y: [], z: [] };
        const t = Vec3();

        const { extraRadius, minRadius, durationMs } = { ...DefaultCameraFocusOptions, ...options };
        const radius = Math.max(sphere.radius + extraRadius, minRadius);

        if (radius <= 1e-3) {
            this.focusSphere(sphere, options);
            return;
        }

        const entityType = StructureProperties.entity.type;

        for (const l of lociArray) {
            if (!StructureElement.Loci.is(l)) continue;
            const extended = StructureElement.Loci.extendToRadius(l, radius);
            StructureElement.Loci.forEachLocation(extended, loc => {
                if (entityType(loc) === 'water') return;

                loc.unit.conformation.position(loc.element, t);
                positions.x.push(t[0]);
                positions.y.push(t[1]);
                positions.z.push(t[2]);
            });
        }

        if (positions.x.length === 0) {
            this.focusSphere(sphere, options);
            return;
        }

        const direction = leastObstructedDirection(positions, {
            origin: sphere.center,
            minDistance: 1e-3,
            sigma: sphere.radius,
        });
        if (!direction) {
            this.focusSphere(sphere, options);
            return;
        }

        Vec3.negate(direction, direction);
        const snapshot = canvas3d.camera.getInvariantFocus(sphere.center, radius, options?.up ?? Vec3.unitY, direction);
        canvas3d.requestCameraReset({ durationMs, snapshot });
    }

    private focusLociBase(loci: Loci | Loci[], options?: Partial<CameraFocusOptions>) {
        const sphere = this.getFocusSphere(loci);
        if (sphere) {
            this.focusSphere(sphere, options);
        }
    }

    focusLoci(loci: Loci | Loci[], options?: Partial<CameraFocusLociOptions>) {
        if (options?.optimizeDirection) {
            this.focusLociOptimized(loci, {
                ...options,
                up: options.optimizeDirectionUp,
            });
        } else {
            this.focusLociBase(loci, options);
        }
    }

    focusSpheres<T>(xs: ReadonlyArray<T>, sphere: (t: T) => Sphere3D | undefined, options?: Partial<CameraFocusOptions> & { principalAxes?: PrincipalAxes, positionToFlip?: Vec3 }) {
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

    focusSphere(sphere: Sphere3D, options?: Partial<CameraFocusOptions> & { principalAxes?: PrincipalAxes, positionToFlip?: Vec3 }) {
        const { canvas3d } = this.plugin;
        if (!canvas3d) return;

        const { extraRadius, minRadius, durationMs } = { ...DefaultCameraFocusOptions, ...options };
        const radius = Math.max(sphere.radius + extraRadius, minRadius);

        if (options?.principalAxes) {
            const snapshot = pcaFocus(this.plugin, radius, options as { principalAxes: PrincipalAxes, positionToFlip?: Vec3 });
            this.plugin.canvas3d?.requestCameraReset({ durationMs, snapshot });
        } else {
            const snapshot = canvas3d.camera.getFocus(sphere.center, radius);
            canvas3d.requestCameraReset({ durationMs, snapshot });
        }
    }

    /** Focus on a set of plugin state object cells (if `options.targets` is non-empty) or on the whole scene (if `options.targets` is empty). */
    focusObject(options: PluginState.SnapshotFocusInfo & { minRadius?: number, durationMs?: number }) {
        if (!this.plugin.canvas3d) return;
        const snapshot = getFocusSnapshot(this.plugin, {
            ...options,
            targets: options.targets?.map(t => ({ ...t, extraRadius: t.extraRadius ?? DefaultCameraFocusOptions.extraRadius })),
            minRadius: options.minRadius ?? DefaultCameraFocusOptions.minRadius,
        });
        this.plugin.canvas3d.requestCameraReset({ snapshot, durationMs: options.durationMs ?? DefaultCameraFocusOptions.durationMs });
    }

    /** Align PCA axes of `structures` (default: all loaded structures) to the screen axes. */
    orientAxes(structures?: Structure[], durationMs?: number) {
        if (!this.plugin.canvas3d) return;
        if (!structures) {
            const structCells = this.plugin.state.data.selectQ(q => q.ofType(PluginStateObject.Molecule.Structure));
            const rootStructCells = structCells.filter(cell => cell.obj && !cell.transform.transformer.definition.isDecorator && !cell.obj.data.parent);
            structures = rootStructCells.map(cell => cell.obj?.data).filter(struct => !!struct) as Structure[];
        }
        const { rotation } = structureLayingTransform(structures);
        const newSnapshot = changeCameraRotation(this.plugin.canvas3d.camera.getSnapshot(), rotation);
        this.setSnapshot(newSnapshot, durationMs);
    }

    /** Align Cartesian axes to the screen axes (X right, Y up). */
    resetAxes(durationMs?: number) {
        if (!this.plugin.canvas3d) return;
        const newSnapshot = changeCameraRotation(this.plugin.canvas3d.camera.getSnapshot(), Mat3.Identity);
        this.setSnapshot(newSnapshot, durationMs);
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