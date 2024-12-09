/**
 * Copyright (c) 2024 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Adam Midlik <midlik@gmail.com>
 */

import { Camera } from '../../../mol-canvas3d/camera';
import { GraphicsRenderObject } from '../../../mol-gl/render-object';
import { Sphere3D } from '../../../mol-math/geometry';
import { BoundaryHelper } from '../../../mol-math/geometry/boundary-helper';
import { Vec3 } from '../../../mol-math/linear-algebra';
import { Loci } from '../../../mol-model/loci';
import { Structure } from '../../../mol-model/structure';
import { PluginContext } from '../../../mol-plugin/context';
import { PluginState } from '../../../mol-plugin/state';
import { StateObject, StateTransform } from '../../../mol-state';
import { PluginStateObject } from '../../objects';


/** Return camera snapshot focused on a plugin state object cell (if `targetRef` is defined)
 * or on the whole scene (if `targetRef` is undefined).
 * If `direction` and `up` are not provided, use current camera orientation. */
export function getFocusSnapshot(plugin: PluginContext, options: PluginState.SnapshotFocusInfo & { minRadius?: number }) {
    if (!plugin.canvas3d) return undefined;
    const targetSpheres = options.targets?.map(target => {
        const bounding = (target.targetRef !== undefined) ? getCellBoundingSphere(plugin, target.targetRef) : getPluginBoundingSphere(plugin);
        if (!bounding) return undefined;
        const radius = target.radius ?? bounding.radius * (target.radiusFactor ?? 1) + (target.extraRadius ?? 0);
        return Sphere3D.create(bounding.center, radius);
    }).filter(sphere => sphere !== undefined);
    const mergedSphere = (targetSpheres && targetSpheres.length > 0) ? boundingSphereOfSpheres(targetSpheres) : getPluginBoundingSphere(plugin);
    return snapshotFromSphereAndDirections(plugin.canvas3d.camera, {
        center: mergedSphere.center,
        radius: Math.max(mergedSphere.radius, options.minRadius ?? 0),
        up: options.up,
        direction: options.direction,
    });
}

const _tmpVec = Vec3();

/** Return camera snapshot for focusing a sphere with given `center` and `radius`,
 * while ensuring given view `direction` (aligns with vector position->target)
 * and `up` (aligns with screen Y axis). */
function snapshotFromSphereAndDirections(camera: Camera, options: { center: Vec3, radius: number, direction?: Vec3, up?: Vec3 }): Partial<Camera.Snapshot> {
    // This might seem to repeat `plugin.canvas3d.camera.getFocus` but avoid flipping
    const target = options.center;
    const radius = Math.max(options.radius, 0.01);
    const direction = options.direction ?? Vec3.sub(Vec3(), camera.target, camera.position);
    const up = Vec3.orthogonalize(Vec3(), direction, options.up ?? camera.up);
    const distance = camera.getTargetDistance(radius);
    const deltaDirection = Vec3.setMagnitude(_tmpVec, direction, distance);
    const position = Vec3.sub(Vec3(), target, deltaDirection);
    return { target, position, up, radius };
}

/** Return the bounding sphere of a plugin state object cell */
export function getCellBoundingSphere(plugin: PluginContext, cellRef: StateTransform.Ref): Sphere3D | undefined {
    const spheres = collectCellBoundingSpheres([], plugin, cellRef);
    if (spheres.length === 0) return undefined;
    if (spheres.length === 1) return spheres[0];
    return boundingSphereOfSpheres(spheres);
}

/** Push bounding spheres within cell `cellRef` to `out`. If a cell does not define bounding spheres, collect bounding spheres from subtree. */
function collectCellBoundingSpheres(out: Sphere3D[], plugin: PluginContext, cellRef: StateTransform.Ref): Sphere3D[] {
    const cell = plugin.state.data.cells.get(cellRef);
    const spheres = getStateObjectBoundingSpheres(cell?.obj);
    if (spheres) {
        out.push(...spheres);
    } else {
        const children = plugin.state.data.tree.children.get(cellRef);
        children.forEach(child => collectCellBoundingSpheres(out, plugin, child));
    }
    return out;
}

/** Return a set of bounding spheres of a plugin state object. Return `undefined` if this plugin state object type does not define bounding spheres. */
function getStateObjectBoundingSpheres(obj: StateObject | undefined): Sphere3D[] | undefined {
    if (!obj) return undefined;
    if (!obj.data) {
        console.warn('Focus: no data');
        return undefined;
    }
    if (obj.data instanceof Structure) {
        const sphere = Loci.getBoundingSphere(Structure.Loci(obj.data));
        return sphere ? [sphere] : [];
    } else if (PluginStateObject.isRepresentation3D(obj)) {
        const out: Sphere3D[] = [];
        for (const renderObject of obj.data.repr.renderObjects) {
            const sphere = renderObject.values.boundingSphere.ref.value;
            if (sphere.radius > 0) out.push(sphere);
        }
        return out;
    }
    return undefined;
}

/** Return the bounding sphere of the whole visible scene. */
export function getPluginBoundingSphere(plugin: PluginContext) {
    const renderObjects = getRenderObjects(plugin, false);
    const spheres = renderObjects.map(r => r.values.boundingSphere.ref.value).filter(sphere => sphere.radius > 0);
    return boundingSphereOfSpheres(spheres);
}

function getRenderObjects(plugin: PluginContext, includeHidden: boolean): GraphicsRenderObject[] {
    let reprCells = Array.from(plugin.state.data.cells.values()).filter(cell => cell.obj && PluginStateObject.isRepresentation3D(cell.obj));
    if (!includeHidden) reprCells = reprCells.filter(cell => !cell.state.isHidden);
    const renderables = reprCells.flatMap(cell => cell.obj!.data.repr.renderObjects);
    return renderables;
}

let boundaryHelper: BoundaryHelper | undefined = undefined;

function boundingSphereOfSpheres(spheres: Sphere3D[]): Sphere3D {
    boundaryHelper ??= new BoundaryHelper('98');
    boundaryHelper.reset();
    for (const s of spheres) boundaryHelper.includeSphere(s);
    boundaryHelper.finishedIncludeStep();
    for (const s of spheres) boundaryHelper.radiusSphere(s);
    return boundaryHelper.getSphere();
}
