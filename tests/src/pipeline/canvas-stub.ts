/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

// The canvas of the plugin benches (many.ts, playback.ts): there is no WebGL in Node, so `plugin.canvas3d` is a stub
// that does what Canvas3D does with representations (add, remove, mark the representations of a loci's root, update) but has no
// scene. It counts the calls that reach the scene, whose cost (gl internals) is measured in the browser.

import { Camera } from '@molstar/graphics/canvas3d/camera';
import { type Canvas3DCameraResetOptions, Canvas3DParams } from '@molstar/graphics/canvas3d/canvas3d';
import { MarkIndex } from '@molstar/graphics/canvas3d/mark-index';
import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import { ValueCell } from '@molstar/core/util/value-cell';
import type { GraphicsRenderObject } from '@molstar/graphics/gl/render-object';
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
 * scene. `sceneUpdated` counts the renderables a real scene updates after a mark that changed something: those whose
 * values changed since the last update (`Scene.updateChanged`, observed as the scene's draw table does). A camera reset with a duration starts a transition; `cameraRestarts` counts the resets that land while the
 * previous transition still runs (by the clock `now`, which playback.ts sets from its frame clock).
 */
export function stubCanvas(graphics: any) {
  const reprs = new Map<Representation.Any, number>();
  /** as Canvas3D's: the representations a mark without one visits */
  const markIndex = new MarkIndex();
  const counts = {
    add: 0,
    remove: 0,
    update: 0,
    mark: 0,
    /** `repr.mark` calls of the marks without a representation (`MarkIndex`) */
    reprMark: 0,
    markChanged: 0,
    sceneUpdated: 0,
    syncVisibility: 0,
    commit: 0,
    setProps: 0,
    cameraReset: 0,
    cameraRestarts: 0,
  };
  let renderObjects = 0;
  /** the representations' render objects, as the scene holds them */
  const objectsOf = new Map<Representation.Any, readonly GraphicsRenderObject[]>();
  const live = new Set<GraphicsRenderObject>();
  /** the objects whose value cells changed since the last update (`Scene.updateChanged` takes them) */
  const changedObjects = new Set<GraphicsRenderObject>();
  const observed = new WeakSet<GraphicsRenderObject>();
  const observe = (ro: GraphicsRenderObject) => {
    if (observed.has(ro)) return;
    observed.add(ro);
    const observer: ValueCell.Observer = {
      get changed() {
        return false;
      },
      set changed(v: boolean) {
        if (v) changedObjects.add(ro);
      },
    };
    for (const cell of Object.values(ro.values as Record<string, any>)) {
      if (cell?.ref) ValueCell.observe(cell, observer);
    }
  };
  const setObjects = (repr: Representation.Any, objects: readonly GraphicsRenderObject[]) => {
    for (const ro of objectsOf.get(repr) ?? []) live.delete(ro);
    if (objects.length === 0) objectsOf.delete(repr);
    else objectsOf.set(repr, objects);
    for (const ro of objects) {
      observe(ro);
      live.add(ro);
    }
  };
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
      // Canvas3D's `add` ends with `scene.update` of the representation's objects: those in the scene not updated again
      for (const ro of repr.renderObjects) if (live.has(ro)) changedObjects.delete(ro);
      setObjects(repr, [...repr.renderObjects]);
      markIndex.set(repr);
    },
    remove(repr: Representation.Any) {
      counts.remove++;
      renderObjects -= reprs.get(repr) ?? 0;
      reprs.delete(repr);
      setObjects(repr, []);
      markIndex.delete(repr);
    },
    update(repr?: Representation.Any) {
      counts.update++;
      // all updated: so are the changed ones; a representation's objects are not updated again (`Scene.update`)
      if (!repr) {
        changedObjects.clear();
        for (const r of reprs.keys()) markIndex.set(r);
      } else {
        for (const ro of repr.renderObjects) changedObjects.delete(ro);
        if (reprs.has(repr)) markIndex.set(repr);
      }
    },
    mark(reprLoci: Representation.Loci, action: MarkerAction) {
      counts.mark++;
      let changed = false;
      if (reprLoci.repr) {
        changed = reprLoci.repr.mark(reprLoci.loci, action);
        if (changed) markIndex.note(reprLoci.repr, action);
      } else
        markIndex.forEach(reprLoci.loci, action, (r) => {
          counts.reprMark++;
          if (!r.mark(reprLoci.loci, action)) return;
          markIndex.note(r, action);
          changed = true;
        });
      if (changed) {
        counts.markChanged++;
        for (const ro of changedObjects) if (live.has(ro)) counts.sceneUpdated++;
        changedObjects.clear();
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
  /** the value cells' versions at the last `changedBytes` call */
  const versions = new WeakMap<object, number>();
  /**
   * Bytes of the typed arrays in the representations' render objects whose value cell changed since the last call:
   * what the renderer would upload again (a proxy: the GPU is not here)
   */
  function changedBytes() {
    let bytes = 0;
    for (const repr of reprs.keys()) {
      for (const ro of repr.renderObjects) {
        const values = ro.values as Record<string, any>;
        for (const k in values) {
          const cell = values[k];
          if (!cell?.ref) continue;
          const v = cell.ref.value;
          const array = ArrayBuffer.isView(v) ? v : v && ArrayBuffer.isView(v.array) ? v.array : undefined;
          if (!array) continue;
          if (versions.get(cell) !== cell.ref.version) bytes += array.byteLength;
          versions.set(cell, cell.ref.version);
        }
      }
    }
    return bytes;
  }
  return {
    canvas,
    counts,
    reprs,
    clock,
    changedBytes,
    get renderObjects() {
      return renderObjects;
    },
  };
}
