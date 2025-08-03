/**
 * Copyright (c) 2023-2024 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Adam Midlik <midlik@gmail.com>
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Camera } from '../../mol-canvas3d/camera';
import { CameraFogParams, Canvas3DParams, Canvas3DProps } from '../../mol-canvas3d/canvas3d';
import { BloomParams } from '../../mol-canvas3d/passes/bloom';
import { DofParams } from '../../mol-canvas3d/passes/dof';
import { OutlineParams } from '../../mol-canvas3d/passes/outline';
import { ShadowParams } from '../../mol-canvas3d/passes/shadow';
import { SsaoParams } from '../../mol-canvas3d/passes/ssao';
import { Vec3 } from '../../mol-math/linear-algebra';
import { getFocusSnapshot, getPluginBoundingSphere } from '../../mol-plugin-state/manager/focus-camera/focus-object';
import { PluginCommands } from '../../mol-plugin/commands';
import { PluginContext } from '../../mol-plugin/context';
import { PluginState } from '../../mol-plugin/state';
import { StateObjectSelector } from '../../mol-state';
import { fovAdjustedPosition } from '../../mol-util/camera';
import { ColorNames } from '../../mol-util/color/names';
import { ParamDefinition } from '../../mol-util/param-definition';
import { decodeColor } from './helpers/utils';
import { MolstarLoadingContext } from './load';
import { SnapshotMetadata } from './mvs-data';
import { MolstarNode, MolstarNodeParams } from './tree/molstar/molstar-tree';
import { MVSTreeSchema } from './tree/mvs/mvs-tree';


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
  * Orient the camera based on a focus node params. **/
export async function setFocus(plugin: PluginContext, focuses: { target: StateObjectSelector, params: MolstarNodeParams<'focus'> }[]) {
    const snapshot = getFocusSnapshot(plugin, {
        ...snapshotFocusInfoFromMvsFocuses(focuses),
        minRadius: DefaultFocusOptions.minRadius,
    });
    if (!snapshot) return;
    resetSceneRadiusFactor(plugin);
    await PluginCommands.Camera.SetSnapshot(plugin, { snapshot });
}

function snapshotFocusInfoFromMvsFocuses(focuses: { target: StateObjectSelector, params: MolstarNodeParams<'focus'> }[]): PluginState.SnapshotFocusInfo {
    const lastFocus = (focuses.length > 0) ? focuses[focuses.length - 1] : undefined;
    const direction = lastFocus?.params.direction ?? MVSTreeSchema.nodes.focus.params.fields.direction.default;
    const up = lastFocus?.params.up ?? MVSTreeSchema.nodes.focus.params.fields.up.default;
    return {
        targets: focuses.map<PluginState.SnapshotFocusTargetInfo>(f => ({
            targetRef: f.target.ref === '-=root=-' ? undefined : f.target.ref, // need to treat root separately so it does not include invisible structure parts etc.
            radius: f.params.radius ?? undefined,
            radiusFactor: f.params.radius_factor,
            extraRadius: f.params.radius_extent,
        })),
        direction: Vec3.create(...direction),
        up: Vec3.create(...up),
    };
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

/** Create object for PluginState.Snapshot.camera based on tree loading context and MVS snapshot metadata */
export function createPluginStateSnapshotCamera(plugin: PluginContext, context: MolstarLoadingContext, metadata: SnapshotMetadata & { previousTransitionDurationMs?: number }): PluginState.Snapshot['camera'] {
    const camera: PluginState.Snapshot['camera'] = {
        transitionStyle: 'animate',
        transitionDurationInMs: metadata.previousTransitionDurationMs ?? 0,
    };
    if (context.camera.cameraParams !== undefined) {
        const currentCameraSnapshot = plugin.canvas3d!.camera.getSnapshot();
        const cameraSnapshot = cameraParamsToCameraSnapshot(plugin, context.camera.cameraParams);
        camera.current = { ...currentCameraSnapshot, ...cameraSnapshot };
    } else {
        camera.focus = snapshotFocusInfoFromMvsFocuses(context.camera.focuses);
    }
    return camera;
}

/** Set canvas properties based on a canvas node. */
export function setCanvas(plugin: PluginContext, node: MolstarNode<'canvas'> | undefined) {
    plugin.canvas3d?.setProps(old => modifyCanvasProps(old, node));
}

/** Create a deep copy of `oldCanvasProps` with values modified according to a canvas node params. */
export function modifyCanvasProps(oldCanvasProps: Canvas3DProps, canvasNode: MolstarNode<'canvas'> | undefined, custom?: Record<string, any>): Canvas3DProps {
    const params = canvasNode?.params;
    const backgroundColor = decodeColor(params?.background_color) ?? DefaultCanvasBackgroundColor;

    const molstar_postprocessing = canvasNode?.custom?.molstar_postprocessing;

    const outline = !!molstar_postprocessing?.enable_outline;
    const outlineParams = molstar_postprocessing?.outline_params;

    const shadow = !!molstar_postprocessing?.enable_shadow;
    const shadowParams = molstar_postprocessing?.shadow_params;

    const occlusion = !!molstar_postprocessing?.enable_ssao;
    const occlusionParams = molstar_postprocessing?.ssao_params;

    const fog = !!molstar_postprocessing?.enable_fog;
    const fogParams = molstar_postprocessing?.fog_params;

    const dof = !!molstar_postprocessing?.enable_depth_of_field;
    const dofParams = molstar_postprocessing?.depth_of_field_params;

    const bloom = !!molstar_postprocessing?.enable_bloom;
    const bloomParams = molstar_postprocessing?.bloom_params;

    return {
        ...oldCanvasProps,
        postprocessing: {
            ...oldCanvasProps.postprocessing,
            outline: outline
                ? { name: 'on', params: { ...ParamDefinition.getDefaultValues(OutlineParams), ...outlineParams } }
                : oldCanvasProps.postprocessing.outline,
            shadow: shadow
                ? { name: 'on', params: { ...ParamDefinition.getDefaultValues(ShadowParams), ...shadowParams } }
                : oldCanvasProps.postprocessing.shadow,
            occlusion: occlusion
                ? { name: 'on', params: { ...ParamDefinition.getDefaultValues(SsaoParams), ...occlusionParams } }
                : oldCanvasProps.postprocessing.occlusion,
            dof: dof
                ? { name: 'on', params: { ...ParamDefinition.getDefaultValues(DofParams), ...dofParams } }
                : oldCanvasProps.postprocessing.dof,
            bloom: bloom
                ? { name: 'on', params: { ...ParamDefinition.getDefaultValues(BloomParams), ...bloomParams } }
                : oldCanvasProps.postprocessing.bloom,
        },
        cameraFog: fog
            ? { name: 'on', params: { ...ParamDefinition.getDefaultValues(CameraFogParams), ...fogParams } }
            : oldCanvasProps.cameraFog,
        renderer: {
            ...oldCanvasProps.renderer,
            backgroundColor: backgroundColor,
        },
    };
}
