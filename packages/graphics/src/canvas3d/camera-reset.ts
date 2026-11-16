/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Sphere3D } from '@molstar/core/math/geometry';
import { Vec3 } from '@molstar/core/math/linear-algebra';
import type { Canvas3DCameraResetOptions } from './canvas3d.js';
import type { TrackballControls } from './controls/trackball.js';
import type { Canvas3DParts } from './parts.js';

/** What the camera reset reads besides the canvas's parts */
export interface Canvas3DCameraResetHost {
  readonly parts: Canvas3DParts;
  readonly controls: TrackballControls;
  /** whether an XR session is presenting (its scale is kept then) */
  isXRPresenting(): boolean;
  /** resets the XR scale (outside a session) */
  resetXRScale(): void;
  /** the scene's radius as the camera's `radiusMax` takes it (`Canvas3DProps.sceneRadiusFactor`) */
  getSceneRadius(): number;
}

/**
 * The camera reset of a canvas: requested by the API (`Canvas3D.requestCameraReset`) or by a scene commit that moved
 * objects out of view, resolved once the whole scene is committed (`resolve`), with the options of the request
 */
export interface Canvas3DCameraReset {
  /** a reset waits for the scene's commit */
  readonly requested: boolean;
  /** requests a reset with `options` (merged into those of an earlier unresolved request) */
  request(options?: Canvas3DCameraResetOptions): void;
  /** requests a reset without a transition (the camera has no state yet) */
  requestImmediate(): void;
  /** before the main scene commits: the bounding sphere its visible objects had */
  beforeCommit(): void;
  /**
   * After the main scene committed: requests a reset where objects moved out of the old bounds or out of view (unless
   * the props reset by hand), immediate where the scene was empty; `empty`: no representation was in the scene
   */
  afterCommit(empty: boolean): void;
  /** Resets the camera to the visible scene where a reset was requested (after the whole scene committed) */
  resolve(): void;
}

export function createCanvas3DCameraReset(host: Canvas3DCameraResetHost): Canvas3DCameraReset {
  const { parts, controls } = host;
  let requested = false;
  const options: Canvas3DCameraResetOptions = {
    durationMs: undefined,
    snapshot: undefined,
    keyframes: undefined,
    easing: undefined,
    trajectory: undefined,
  };
  /** the bounding sphere of the visible objects before the last commit */
  const oldBoundingSphereVisible = Sphere3D();
  const cameraSphere = Sphere3D();

  function shouldReset() {
    const { camera, scene } = parts;
    if (camera.state.radiusMax === 0) return true;

    if (camera.transition.inTransition || options.snapshot) return false;

    let cameraSphereOverlapsNone = true,
      isEmpty = true;
    Sphere3D.set(cameraSphere, camera.state.target, camera.state.radius);

    // check if any renderable has moved outside of the old bounding sphere
    // and if no renderable is overlapping with the camera sphere
    for (const r of scene.renderables) {
      if (!r.state.visible) continue;

      const b = r.values.boundingSphere.ref.value;
      if (!b.radius) continue;

      isEmpty = false;
      const cameraDist = Vec3.distance(cameraSphere.center, b.center);
      if (
        (cameraDist > cameraSphere.radius || cameraDist > b.radius || b.radius > camera.state.radiusMax) &&
        !Sphere3D.includes(oldBoundingSphereVisible, b)
      )
        return true;
      if (Sphere3D.overlaps(cameraSphere, b)) cameraSphereOverlapsNone = false;
    }

    return cameraSphereOverlapsNone || (!isEmpty && cameraSphere.radius <= 0.1);
  }

  return {
    get requested() {
      return requested;
    },
    request: (next) => {
      Object.assign(options, next);
      requested = true;
    },
    requestImmediate: () => {
      requested = true;
      options.durationMs = 0;
    },
    beforeCommit: () => {
      Sphere3D.copy(oldBoundingSphereVisible, parts.scene.boundingSphereVisible);
    },
    afterCommit: (empty) => {
      if (!parts.props.camera.manualReset && (empty || shouldReset())) requested = true;
      if (oldBoundingSphereVisible.radius === 0) options.durationMs = 0;
    },
    resolve: () => {
      if (!requested) return;
      const { camera, scene, props: p } = parts;

      if (!host.isXRPresenting()) host.resetXRScale();

      const { center, radius } = scene.boundingSphereVisible;

      if (radius > 0) {
        const duration = options.durationMs ?? p.cameraResetDurationMs;
        const easing = options.easing ?? p.cameraResetEasing;
        const trajectory = options.trajectory ?? p.cameraResetTrajectory;
        const autoAdjustControls = controls.props.autoAdjustMinMaxDistance;
        if (autoAdjustControls.name === 'on') {
          const minDistance =
            autoAdjustControls.params.minDistanceFactor * radius + autoAdjustControls.params.minDistancePadding;
          const maxDistance = Math.max(
            autoAdjustControls.params.maxDistanceFactor * radius,
            autoAdjustControls.params.maxDistanceMin,
          );
          controls.setProps({ minDistance, maxDistance });
        }
        const focus = camera.getFocus(center, radius);
        const next = typeof options.snapshot === 'function' ? options.snapshot(scene, camera) : options.snapshot;
        const snapshot = next ? { ...focus, ...next } : focus;
        camera.setState({ ...snapshot, radiusMax: host.getSceneRadius() }, duration, {
          keyframes: options.keyframes,
          easing,
          trajectory,
        });
      }

      options.durationMs = void 0;
      options.snapshot = void 0;
      options.keyframes = void 0;
      options.easing = void 0;
      options.trajectory = void 0;

      requested = false;
    },
  };
}
