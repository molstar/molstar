/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The canvas of the plugin benches (many.ts, playback.ts): there is no WebGL in Node, so `plugin.canvas3d` is a stub
// that does what Canvas3D does with representations (add, remove, mark every representation, update) but has no
// scene. It counts the calls that reach the scene, whose cost (frozen gl internals) is measured in the browser.

import { Camera } from '@molstar/graphics/canvas3d/camera';
import { type Canvas3DCameraResetOptions, Canvas3DParams } from '@molstar/graphics/canvas3d/canvas3d';
import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import type { Representation } from '@molstar/graphics/repr/representation';
import type { MarkerAction } from '@molstar/core/util/marker-action';

/** the parts of a GraphicsContext the CPU geometry paths read (bench.ts): impostors on, no GPU compute */
export function graphicsStub(meshMode: boolean) {
  return meshMode
    ? undefined
    : ({
        capabilities: { fragDepth: true, flatInterpolation: true, texture: true },
        resources: {
          limits: { maxTextureSize: 16384, max3dTextureSize: 2048 },
          texture: () => new Proxy({} as any, { get: (t, k) => (k === 'then' ? undefined : k in t ? t[k] : () => {}) }),
          ownsTexture: () => true,
        },
        namedTextures: {},
      } as any);
}

/**
 * What Canvas3D does with representations (canvas3d.ts `add`, `remove`, `mark`/`resolveMarking`, `update`), without a
 * scene. `scene*` count the renderables a real scene would walk: `scene.update(void 0)` after a mark goes over all of
 * them. A camera reset with a duration starts a transition; `cameraRestarts` counts the resets that land while the
 * previous transition still runs (by the clock `now`, which playback.ts sets from its frame clock).
 */
export function stubCanvas(graphics: any) {
  const reprs = new Map<Representation.Any, number>();
  const counts = {
    add: 0,
    remove: 0,
    update: 0,
    mark: 0,
    markChanged: 0,
    sceneUpdateAll: 0,
    syncVisibility: 0,
    commit: 0,
    setProps: 0,
    cameraReset: 0,
    cameraRestarts: 0,
  };
  let renderObjects = 0;
  let transitionEnd = -Infinity;
  const clock = { now: 0 };
  const camera = new Camera();
  const canvas = {
    graphics,
    camera,
    props: PD.getDefaultValues(Canvas3DParams),
    add(repr: Representation.Any) {
      counts.add++;
      renderObjects += repr.renderObjects.length - (reprs.get(repr) ?? 0);
      reprs.set(repr, repr.renderObjects.length);
    },
    remove(repr: Representation.Any) {
      counts.remove++;
      renderObjects -= reprs.get(repr) ?? 0;
      reprs.delete(repr);
    },
    update() {
      counts.update++;
    },
    mark(reprLoci: Representation.Loci, action: MarkerAction) {
      counts.mark++;
      let changed = false;
      if (reprLoci.repr) changed = reprLoci.repr.mark(reprLoci.loci, action);
      else for (const r of reprs.keys()) changed = r.mark(reprLoci.loci, action) || changed;
      if (changed) {
        counts.markChanged++;
        counts.sceneUpdateAll += renderObjects;
      }
    },
    syncVisibility() {
      counts.syncVisibility++;
    },
    commit() {
      counts.commit++;
    },
    requestDraw() {},
    /** shader programs are compiled by gl: nothing to do here */
    async precompile() {},
    dispose() {},
    requestCameraReset(options?: Canvas3DCameraResetOptions) {
      counts.cameraReset++;
      const duration = options?.durationMs ?? 0;
      if (duration > 0 && clock.now < transitionEnd) counts.cameraRestarts++;
      transitionEnd = duration > 0 ? clock.now + duration : -Infinity;
      if (options?.snapshot) camera.setState(options.snapshot as Camera.Snapshot, 0);
    },
    setProps() {
      counts.setProps++;
    },
    get boundingSphere() {
      return { center: [0, 0, 0], radius: 0 };
    },
  };
  return {
    canvas,
    counts,
    reprs,
    clock,
    get renderObjects() {
      return renderObjects;
    },
  };
}
