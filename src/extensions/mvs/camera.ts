/**
 * Copyright (c) 2023-2024 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Adam Midlik <midlik@gmail.com>
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Camera } from '../../mol-canvas3d/camera';
import { Canvas3DParams } from '../../mol-canvas3d/canvas3d';
import { Vec3 } from '../../mol-math/linear-algebra';
import { getFocusSnapshot, getPluginBoundingSphere } from '../../mol-plugin-state/manager/focus-camera/focus-object';
import { PluginCommands } from '../../mol-plugin/commands';
import { PluginContext } from '../../mol-plugin/context';
import { StateObjectSelector } from '../../mol-state';
import { ColorNames } from '../../mol-util/color/names';
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
    const radius = Vec3.distance(target, position) / 2;
    if (plugin.canvas3d) position = fovAdjustedPosition(target, position, plugin.canvas3d.camera.state.mode, plugin.canvas3d.camera.state.fov);
    const up = Vec3.create(...params.up);
    Vec3.orthogonalize(up, Vec3.sub(_tmpVec, target, position), up);
    const snapshot: Partial<Camera.Snapshot> = { target, position, up, radius, radiusMax: radius };
    return snapshot;
}

/** Focus the camera on the bounding sphere of a (sub)structure (or on the whole scene if `structureNodeSelector` is undefined).
  * Orient the camera based on a focus node params.
  **/
export async function setFocus(plugin: PluginContext, structureNodeSelector: StateObjectSelector | undefined, params: MolstarNodeParams<'focus'> = MVSDefaults.focus) {
    const snapshot = getFocusSnapshot(plugin, structureNodeSelector?.ref, {
        direction: Vec3.create(...params.direction),
        up: Vec3.create(...params.up),
        extraRadius: DefaultFocusOptions.extraRadius,
        minRadius: DefaultFocusOptions.minRadius,
    });
    if (!snapshot) return;
    resetSceneRadiusFactor(plugin);
    await PluginCommands.Camera.SetSnapshot(plugin, { snapshot });
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
