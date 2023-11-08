/**
 * Copyright (c) 2023 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Adam Midlik <midlik@gmail.com>
 */

import { Camera } from '../../mol-canvas3d/camera';
import { GraphicsRenderObject } from '../../mol-gl/render-object';
import { Sphere3D } from '../../mol-math/geometry';
import { BoundaryHelper } from '../../mol-math/geometry/boundary-helper';
import { Vec3 } from '../../mol-math/linear-algebra';
import { Loci } from '../../mol-model/loci';
import { Structure } from '../../mol-model/structure';
import { PluginStateObject } from '../../mol-plugin-state/objects';
import { PluginCommands } from '../../mol-plugin/commands';
import { PluginContext } from '../../mol-plugin/context';
import { StateObjectSelector } from '../../mol-state';
import { ColorNames } from '../../mol-util/color/names';

import { decodeColor } from './helpers/utils';
import { ParamsOfKind } from './tree/generic/tree-schema';
import { MolstarTree } from './tree/molstar/molstar-tree';
import { MVSDefaults } from './tree/mvs/mvs-defaults';


/** Defined in `../../mol-plugin-state/manager/camera.ts` but private */
const DefaultCameraFocusOptions = {
    minRadius: 5,
    extraRadiusForFocus: 4,
    extraRadiusForZoomAll: 0,
};
const DefaultCanvasBackgroundColor = ColorNames.white;


/** Set the camera based on a camera node params. */
export async function setCamera(plugin: PluginContext, params: ParamsOfKind<MolstarTree, 'camera'>) {
    const target = Vec3.create(...params.target);
    const position = Vec3.create(...params.position);
    const up = Vec3.create(...params.up);
    const snapshot: Partial<Camera.Snapshot> = { target, position, up };
    await PluginCommands.Camera.SetSnapshot(plugin, { snapshot });
}

/** Focus the camera on the bounding sphere of a (sub)structure (or on the whole scene if `structureNodeSelector` is null).
 * Orient the camera based on a focus node params. */
export async function setFocus(plugin: PluginContext, structureNodeSelector: StateObjectSelector | undefined, params: ParamsOfKind<MolstarTree, 'focus'> = MVSDefaults.focus) {
    let structure: Structure | undefined = undefined;
    if (structureNodeSelector) {
        const cell = plugin.state.data.cells.get(structureNodeSelector.ref);
        structure = cell?.obj?.data;
        if (!structure) console.warn('Focus: no structure');
        if (!(structure instanceof Structure)) {
            console.warn('Focus: cannot apply to a non-structure node');
            structure = undefined;
        }
    }
    const boundingSphere = structure ? Loci.getBoundingSphere(Structure.Loci(structure)) : getPluginBoundingSphere(plugin);
    if (boundingSphere && plugin.canvas3d) {
        // cannot use plugin.canvas3d.camera.getFocus with up+direction, because it sometimes flips orientation
        // await PluginCommands.Camera.Focus(plugin, { center: boundingSphere.center, radius: boundingSphere.radius }); // this could not set orientation
        const target = boundingSphere.center;
        const extraRadius = structure ? DefaultCameraFocusOptions.extraRadiusForFocus : DefaultCameraFocusOptions.extraRadiusForZoomAll;
        const sphereRadius = Math.max(boundingSphere.radius + extraRadius, DefaultCameraFocusOptions.minRadius);
        const distance = getFocusDistance(plugin.canvas3d.camera, boundingSphere.center, sphereRadius) ?? 100;
        const direction = Vec3.create(...params.direction);
        Vec3.setMagnitude(direction, direction, distance);
        const position = Vec3.sub(Vec3(), target, direction);
        const up = Vec3.create(...params.up);
        const snapshot: Partial<Camera.Snapshot> = { target, position, up, radius: sphereRadius };
        await PluginCommands.Camera.SetSnapshot(plugin, { snapshot });
    }
}

/** Calculate the necessary distance between the camera position and center of the target,
 * if we want to zoom the target with the given radius. */
function getFocusDistance(camera: Camera, target: Vec3, radius: number) {
    const p = camera.getFocus(target, radius);
    if (!p.position || !p.target) return undefined;
    return Vec3.distance(p.position, p.target);
}

/** Compute the bounding sphere of the whole scene. */
function getPluginBoundingSphere(plugin: PluginContext) {
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

/** Set canvas properties based on a canvas node params. */
export function setCanvas(plugin: PluginContext, params: ParamsOfKind<MolstarTree, 'canvas'> | undefined) {
    const backgroundColor = decodeColor(params?.background_color) ?? DefaultCanvasBackgroundColor;
    if (backgroundColor !== plugin.canvas3d?.props.renderer.backgroundColor) {
        plugin.canvas3d?.setProps(old => ({
            ...old,
            renderer: {
                ...old.renderer,
                backgroundColor: backgroundColor,
            }
        }));
    }
}
