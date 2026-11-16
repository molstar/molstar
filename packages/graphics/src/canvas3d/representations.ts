/**
 * Copyright (c) 2018-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

// A canvas's representations: the render objects each has in the scene of its layer (the scene, the overlay, a screen
// scene), kept in step as it updates (`Representation.updated`), moved with its layer, and dropped by a context change
// until it updated them

import type { PrecompileRecorder } from '@molstar/graphics/gl/precompile';
import { type GraphicsRenderObject, hasForeignResources } from '@molstar/graphics/gl/render-object';
import type { GraphicsRenderable } from '@molstar/graphics/gl/renderable';
import type { Scene } from '@molstar/graphics/gl/scene';
import type { Representation } from '@molstar/graphics/repr/representation';
import { ValueCell } from '@molstar/core/util';
import { SetUtils } from '@molstar/core/util/set';
import type { Subscription } from 'rxjs';
import { MarkIndex } from './mark-index.js';
import type { Canvas3DParts } from './parts.js';

/** What the representations call back into the canvas */
export interface Canvas3DRepresentationsHost {
  readonly parts: Canvas3DParts;
  readonly precompileRecorder: PrecompileRecorder;
  /** the canvas has work (wakes its loop) */
  busy: () => void;
  /** frames wait for the pipelines of renderables updated in place (`Canvas3D.holdFramesUntilPrepared`) */
  holdFramesUntilPrepared: (inFlight: boolean) => void;
  /** `owner`'s objects were added or removed (`Canvas3DPreparer.markPreparedObjectsAdded`) */
  markPreparedObjectsAdded: (owner: object) => void;
  /** what was added, removed or updated is drawn once all is committed */
  changed: () => void;
  /** debug mode: the stats logged after a change */
  logStats: () => void;
}

export interface Canvas3DRepresentations {
  /** the representations a loci marked on all of them may change (`applyMark`) */
  readonly markIndex: MarkIndex;
  /** Adds `repr`'s render objects to the scene of its layer, or brings them up to date (`Canvas3D.add`) */
  add: (repr: Representation.Any) => void;
  /** Takes `repr`'s render objects out (`Canvas3D.remove`) */
  remove: (repr: Representation.Any) => void;
  /**
   * The values of `repr`'s render objects updated in place, all representations' without one (`Canvas3D.update`);
   * false where nothing was updated in place (`repr` not added, or added again for its new layer)
   */
  update: (repr: Representation.Any | undefined, keepSphere: boolean) => boolean;
  /** Removes all representations and the objects of the scene, the overlay and the debug helper (`Canvas3D.clear`) */
  clear: () => void;
  has: (repr: Representation.Any) => boolean;
  /** the representations added */
  keys: () => IterableIterator<Representation.Any>;
  /** the scene of the layer `repr`'s render objects are in, undefined where not added */
  sceneOf: (repr: Representation.Any) => Scene | undefined;
  /** the representations of the scene layer: what the scene shows (gizmos and widgets of the other layers not) */
  getSceneCount: () => number;
  /** whether a representation added holds `o` */
  owns: (o: GraphicsRenderObject) => boolean;
  /** `o` left out of its scene by a context change, until its representation updates it */
  drop: (o: GraphicsRenderObject) => void;
  isDropped: (o: GraphicsRenderObject) => boolean;
}

export function createCanvas3DRepresentations(host: Canvas3DRepresentationsHost): Canvas3DRepresentations {
  const { parts, precompileRecorder } = host;
  const reprRenderObjects = new Map<Representation.Any, Set<GraphicsRenderObject>>();
  /** the scene of the layer each representation's render objects are in, see `layerScene` */
  const reprScenes = new Map<Representation.Any, Scene>();
  const reprUpdatedSubscriptions = new Map<Representation.Any, Subscription>();
  const markIndex = new MarkIndex();
  /** render objects left out of the scene by a context change until their representation updated them */
  const droppedRenderObjects = new Set<GraphicsRenderObject>();

  /** The scene of the layer `repr` is drawn in (`Representation.State.layer`) */
  function layerScene(repr: Representation.Any): Scene {
    const { helper } = parts;
    const { layer } = repr.state;
    if (layer === 'overlay') return helper.overlay.scene;
    if (layer === 'screen') return helper.screen.getScene(repr, repr.state.screen);
    return parts.scene;
  }

  /** Takes the render objects of `repr` out of their layer's scene (a screen scene of its own goes with them) */
  function removeFromLayer(repr: Representation.Any, target: Scene, renderObjects: Set<GraphicsRenderObject>) {
    const { helper } = parts;
    if (helper.screen.has(target)) {
      helper.screen.release(repr);
      renderObjects.forEach((o) => droppedRenderObjects.delete(o));
    } else {
      removeRenderObjects(target, renderObjects);
    }
  }

  function removeRenderObjects(target: Scene, renderObjects: Set<GraphicsRenderObject>) {
    renderObjects.forEach((o) => {
      target.remove(o);
      droppedRenderObjects.delete(o);
    });
  }

  /**
   * The representation's state that its render objects' values carry (`Representation.State.occludedAlpha`,
   * `pickAtAnyAlpha`): on what it was made with, again after a `setState` (`update`) or a rebuild (`add`)
   */
  function applyStateValues(repr: Representation.Any) {
    const { occludedAlpha, pickAtAnyAlpha } = repr.state;
    const occluded = Math.min(Math.max(occludedAlpha, 0), 1);
    const pickAlphaMin = pickAtAnyAlpha ? 1 : 0;
    const { renderObjects } = repr;
    for (let i = 0, il = renderObjects.length; i < il; ++i) {
      const { values } = renderObjects[i];
      ValueCell.updateIfChanged(values.uOccludedAlpha, occluded);
      ValueCell.updateIfChanged(values.uPickingAlphaMin, pickAlphaMin);
    }
  }

  function add(repr: Representation.Any) {
    registerAutoUpdate(repr);
    host.busy();
    host.markPreparedObjectsAdded(repr);

    const { scene, graphics } = parts;
    const target = layerScene(repr);
    let oldRO = reprRenderObjects.get(repr);
    const oldScene = reprScenes.get(repr);
    if (oldRO && oldScene && oldScene !== target) {
      // moved to another layer
      removeFromLayer(repr, oldScene, oldRO);
      oldRO = undefined;
    }
    const newRO = new Set<GraphicsRenderObject>();
    const overlay = target !== scene;
    repr.renderObjects.forEach((o) => {
      newRO.add(o);
      // the overlay's occlusion (`apply_occlusion`), before its renderable is created with it
      ValueCell.updateIfChanged(o.values.dOverlay, overlay);
    });
    applyStateValues(repr);

    if (oldRO) {
      const old = oldRO;
      if (!SetUtils.areEqual(newRO, old)) {
        newRO.forEach((o) => {
          if (!old.has(o)) target.add(o);
        });
        old.forEach((o) => {
          if (!newRO.has(o)) target.remove(o);
        });
      }
      if (droppedRenderObjects.size > 0) {
        old.forEach((o) => {
          if (!droppedRenderObjects.has(o)) return;
          if (!newRO.has(o)) {
            droppedRenderObjects.delete(o);
          } else if (!hasForeignResources(o, graphics)) {
            droppedRenderObjects.delete(o);
            target.add(o);
          }
        });
      }
    } else {
      repr.renderObjects.forEach((o) => target.add(o));
    }
    reprRenderObjects.set(repr, newRO);
    reprScenes.set(repr, target);
    markIndex.set(repr);

    target.update(repr.renderObjects, false);
    // the objects kept were updated in place (e.g. a theme change)
    if (oldRO) holdForUpdated(repr.renderObjects, target);
    host.changed();
    host.logStats();
  }

  /**
   * Records the committed renderables of `renderObjects` again after their values were updated in place; where their
   * new defines (e.g. another color type) may need programs not compiled yet, frames wait for them as after a props
   * change, instead of skipping their draws meanwhile (WebGPU)
   */
  function holdForUpdated(renderObjects: ArrayLike<GraphicsRenderObject>, target: Scene) {
    const renderables: GraphicsRenderable[] = [];
    for (let i = 0, il = renderObjects.length; i < il; ++i) {
      const renderable = target.get(renderObjects[i]);
      if (renderable) renderables.push(renderable);
    }
    if (precompileRecorder.recordAgain(renderables)) host.holdFramesUntilPrepared(true);
  }

  function remove(repr: Representation.Any) {
    unregisterAutoUpdate(repr);
    host.busy();
    host.markPreparedObjectsAdded(repr);
    markIndex.delete(repr);

    const renderObjects = reprRenderObjects.get(repr);
    if (renderObjects) {
      removeFromLayer(repr, reprScenes.get(repr) ?? parts.scene, renderObjects);
      reprRenderObjects.delete(repr);
      reprScenes.delete(repr);
      host.changed();
      host.logStats();
    }
  }

  function update(repr: Representation.Any | undefined, keepSphere: boolean) {
    // updated values may have changed defines: the updated renderables recorded again
    if (repr) {
      if (!reprRenderObjects.has(repr)) return false;
      if (reprScenes.get(repr) !== layerScene(repr)) {
        // its layer changed (`setState`)
        add(repr);
        return false;
      }
      markIndex.set(repr);
      applyStateValues(repr);
      const { renderObjects } = repr;
      const target = layerScene(repr);
      target.update(renderObjects, keepSphere);
      // only the representation's, not the whole scene's: callers update one representation at a time (a drag)
      holdForUpdated(renderObjects, target);
    } else {
      const { scene, helper } = parts;
      reprRenderObjects.forEach((_, r) => markIndex.set(r));
      scene.update(void 0, keepSphere);
      helper.overlay.scene.update(void 0, keepSphere);
      for (const s of helper.screen.scenes) s.update(void 0, keepSphere);
      precompileRecorder.invalidate();
      let changed = precompileRecorder.recordAgain(scene.renderables);
      if (precompileRecorder.recordAgain(helper.overlay.scene.renderables)) changed = true;
      for (const s of helper.screen.scenes) if (precompileRecorder.recordAgain(s.renderables)) changed = true;
      if (changed) host.holdFramesUntilPrepared(true);
    }
    return true;
  }

  function clear() {
    const { scene, helper } = parts;
    reprUpdatedSubscriptions.forEach((v) => v.unsubscribe());
    reprUpdatedSubscriptions.clear();
    reprRenderObjects.clear();
    // the representations' screen scenes, not the corner axes'
    for (const [repr, s] of reprScenes) if (helper.screen.has(s)) helper.screen.release(repr);
    reprScenes.clear();
    markIndex.clear();
    scene.clear();
    helper.overlay.scene.clear();
    helper.debug.clear();
  }

  function registerAutoUpdate(repr: Representation.Any) {
    if (reprUpdatedSubscriptions.has(repr)) return;

    reprUpdatedSubscriptions.set(
      repr,
      repr.updated.subscribe((_) => {
        if (!repr.state.syncManually) add(repr);
      }),
    );
  }

  function unregisterAutoUpdate(repr: Representation.Any) {
    const updatedSubscription = reprUpdatedSubscriptions.get(repr);
    if (updatedSubscription) {
      updatedSubscription.unsubscribe();
      reprUpdatedSubscriptions.delete(repr);
    }
  }

  return {
    markIndex,
    add,
    remove,
    update,
    clear,
    has: (repr) => reprRenderObjects.has(repr),
    keys: () => reprRenderObjects.keys(),
    sceneOf: (repr) => reprScenes.get(repr),
    getSceneCount: () => {
      let count = 0;
      for (const s of reprScenes.values()) if (s === parts.scene) count += 1;
      return count;
    },
    owns: (o) => {
      for (const renderObjects of reprRenderObjects.values()) {
        if (renderObjects.has(o)) return true;
      }
      return false;
    },
    drop: (o) => {
      droppedRenderObjects.add(o);
    },
    isDropped: (o) => droppedRenderObjects.has(o),
  };
}
