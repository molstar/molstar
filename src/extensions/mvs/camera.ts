/**
 * Copyright (c) 2023-2024 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Adam Midlik <midlik@gmail.com>
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Camera } from '../../mol-canvas3d/camera';
import { Canvas3DParams } from '../../mol-canvas3d/canvas3d';
import { GraphicsRenderObject } from '../../mol-gl/render-object';
import { Sphere3D } from '../../mol-math/geometry';
import { BoundaryHelper } from '../../mol-math/geometry/boundary-helper';
import { Vec3 } from '../../mol-math/linear-algebra';
import { Loci } from '../../mol-model/loci';
import { Structure } from '../../mol-model/structure';
import { PluginStateObject } from '../../mol-plugin-state/objects';
import { PluginCommands } from '../../mol-plugin/commands';
import { PluginContext } from '../../mol-plugin/context';
import { StateObjectSelector, StateTransform } from '../../mol-state';
import { ColorNames } from '../../mol-util/color/names';
import { MVSPrimitivesData } from './components/primitives';
import { decodeColor } from './helpers/utils';
import { MolstarNodeParams } from './tree/molstar/molstar-tree';
import { MVSDefaults } from './tree/mvs/mvs-defaults';


const DefaultFocusOptions = {
    minRadius: 5,
    extraRadius: 0,
};
const DefaultCanvasBackgroundColor = ColorNames.white;


const _tmpVec = Vec3();

/** Set the camera position to the current position (thus suppress automatic adjustment). */
export async function suppressCameraAutoreset(plugin: PluginContext) {
    const snapshot: Partial<Camera.Snapshot> = { ...plugin.canvas3d?.camera.state, radius: Infinity }; // `radius: Infinity` avoids clipping when the scene expands
    adjustSceneRadiusFactor(plugin, snapshot.target);
    await PluginCommands.Camera.SetSnapshot(plugin, { snapshot });
}

/** Set the camera based on a camera node params. */
export async function setCamera(plugin: PluginContext, params: MolstarNodeParams<'camera'>) {
    const snapshot = cameraParamsToCameraSnapshot(plugin, params);
    adjustSceneRadiusFactor(plugin, snapshot.target);
    await PluginCommands.Camera.SetSnapshot(plugin, { snapshot });
}

export function cameraParamsToCameraSnapshot(plugin: PluginContext, params: MolstarNodeParams<'camera'>): Partial<Camera.Snapshot> {
    const target = Vec3.create(...params.target);
    let position = Vec3.create(...params.position);
    if (plugin.canvas3d) position = fovAdjustedPosition(target, position, plugin.canvas3d.camera.state.mode, plugin.canvas3d.camera.state.fov);
    const up = Vec3.create(...params.up);
    Vec3.orthogonalize(up, Vec3.sub(_tmpVec, target, position), up);
    const snapshot: Partial<Camera.Snapshot> = { target, position, up, radius: Infinity }; // `radius: Infinity` avoids clipping (ensures covering the whole scene)
    return snapshot;
}

function getRenderObjectsBoundary(objects: ReadonlyArray<GraphicsRenderObject>) {
    const spheres: Sphere3D[] = [];
    for (const o of objects) {
        const s = o.values.boundingSphere.ref.value;
        if (s.radius === 0) continue;
        spheres.push(s);
    }
    if (spheres.length === 0) return;
    if (spheres.length === 1) return spheres[0];
    return boundingSphereOfSpheres(spheres);
}

/** Focus the camera on the bounding sphere of a (sub)structure (or on the whole scene if `structureNodeSelector` is null).
  * Orient the camera based on a focus node params.
  **/
export async function setFocus(plugin: PluginContext, structureNodeSelector: StateObjectSelector | undefined, params: MolstarNodeParams<'focus'> = MVSDefaults.focus) {
    const snapshot = getFocusSnapshot(plugin, structureNodeSelector?.ref, { direction: Vec3.create(...params.direction), up: Vec3.create(...params.up), extraRadius: DefaultFocusOptions.extraRadius, minRadius: DefaultFocusOptions.minRadius });
    if (!snapshot) return;
    resetSceneRadiusFactor(plugin);
    await PluginCommands.Camera.SetSnapshot(plugin, { snapshot });
}

export function getFocusSnapshot(plugin: PluginContext, targetNodeRef: StateTransform.Ref | undefined, options: { direction?: Vec3, up?: Vec3, extraRadius?: number, minRadius?: number }) {
    if (!plugin.canvas3d) return undefined;
    const boundingSphere = getFocusBoundingSphere(plugin, targetNodeRef);
    if (!boundingSphere) return undefined;
    return snapshotFromSphereAndDirections(plugin.canvas3d.camera, {
        center: boundingSphere.center,
        radius: Math.max(boundingSphere.radius + (options.extraRadius ?? 0), options.minRadius ?? 0),
        up: options.up,
        direction: options.direction,
    });
}

function getFocusBoundingSphere(plugin: PluginContext, targetNodeRef: StateTransform.Ref | undefined) {
    if (targetNodeRef) {
        const cell = plugin.state.data.cells.get(targetNodeRef);
        const data = cell?.obj?.data;
        if (!data) console.warn('Focus: no structure');
        if (data instanceof Structure) {
            return Loci.getBoundingSphere(Structure.Loci(data));
        } else if (PluginStateObject.isRepresentation3D(cell?.obj)) {
            return getRenderObjectsBoundary(cell.obj.data.repr.renderObjects);
        } else if (MVSPrimitivesData.is(cell?.obj)) {
            const representations = plugin.state.data.selectQ(q =>
                q.byRef(cell.transform.ref).subtree().filter(c => PluginStateObject.isRepresentation3D(c?.obj))
            );
            const renderObjects = representations.flatMap(r => r.obj?.data?.repr?.renderObjects ?? []);
            if (renderObjects.length) {
                return getRenderObjectsBoundary(renderObjects);
            }
        } else {
            console.warn('Focus: cannot apply to the specified node type');
        }
    }
    return getPluginBoundingSphere(plugin);
}

/** Adjust `sceneRadiusFactor` property so that the current scene is not cropped */
function adjustSceneRadiusFactor(plugin: PluginContext, cameraTarget: Vec3 | undefined) {
    if (!cameraTarget) return;
    const boundingSphere = getPluginBoundingSphere(plugin);
    const offset = Vec3.distance(cameraTarget, boundingSphere.center);
    const sceneRadiusFactor = boundingSphere.radius > 0 ? ((boundingSphere.radius + offset) / boundingSphere.radius) : 1;
    plugin.canvas3d?.setProps({ sceneRadiusFactor });
}

/** Reset `sceneRadiusFactor` property to the default value */
function resetSceneRadiusFactor(plugin: PluginContext) {
    const sceneRadiusFactor = Canvas3DParams.sceneRadiusFactor.defaultValue;
    plugin.canvas3d?.setProps({ sceneRadiusFactor });
}

/** Return camera snapshot for focusing a sphere with given `center` and `radius`,
 * while ensuring given view `direction` (aligns with vector position->target)
 * and `up` (aligns with screen Y axis). */
function snapshotFromSphereAndDirections(camera: Camera, options: { center: Vec3, radius: number, direction?: Vec3, up?: Vec3 }): Partial<Camera.Snapshot> {
    // This might seem to repeat `plugin.canvas3d.camera.getFocus` but avoid flipping
    const target = options.center;
    const radius = options.radius;
    const direction = options.direction ?? Vec3.sub(Vec3(), camera.target, camera.position);
    const up = Vec3.orthogonalize(Vec3(), direction, options.up ?? camera.up);
    const distance = camera.getTargetDistance(radius);
    const deltaDirection = Vec3.setMagnitude(_tmpVec, direction, distance);
    const position = Vec3.sub(Vec3(), target, deltaDirection);
    return { target, position, up, radius };
}

/** Return the distance adjustment ratio for conversion from the "reference camera"
 * to a camera with an arbitrary field of view `fov`. */
function distanceAdjustment(mode: Camera.Mode, fov: number) {
    if (mode === 'orthographic') return 1 / (2 * Math.tan(fov / 2));
    else return 1 / (2 * Math.sin(fov / 2));
}

/** Return the position for a camera with an arbitrary field of view `fov`
 * necessary to just fit into view the same sphere (with center at `target`)
 * as the "reference camera" placed at `refPosition` would fit, while keeping the camera orientation.
 * The "reference camera" is a camera which can just fit into view a sphere of radius R with center at distance 2R
 * (this corresponds to FOV = 2 * asin(1/2) in perspective mode or FOV = 2 * atan(1/2) in orthographic mode). */
function fovAdjustedPosition(target: Vec3, refPosition: Vec3, mode: Camera.Mode, fov: number) {
    const delta = Vec3.sub(Vec3(), refPosition, target);
    const adjustment = distanceAdjustment(mode, fov);
    return Vec3.scaleAndAdd(delta, target, delta, adjustment); // return target + delta * adjustment
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
export function setCanvas(plugin: PluginContext, params: MolstarNodeParams<'canvas'> | undefined) {
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
