/**
 * Copyright (c) 2018-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { GraphicsContext } from './context.js';
import { type GraphicsRenderObject, createRenderable, transferTextures } from './render-object.js';
import { Object3D } from './object3d.js';
import { Sphere3D } from '@molstar/core/math/geometry/primitives/sphere3d';
import { CommitQueue } from './commit-queue.js';
import { now } from '@molstar/core/util/now';
import { arraySetRemove } from '@molstar/core/util/array';
import { BoundaryHelper } from '@molstar/core/math/geometry/boundary-helper';
import type { GraphicsRenderable } from './renderable.js';
import type { Transparency } from './render-item.js';
import type { GlobalDefines } from './renderable/schema.js';
import { ValueCell } from '@molstar/core/util/value-cell';
import { createDrawTable, DrawFact, DrawRow, DrawState, type DrawTable } from './draw-table.js';

const boundaryHelper = new BoundaryHelper('98');

function calculateBoundingSphere(
  renderables: GraphicsRenderable[],
  boundingSphere: Sphere3D,
  onlyVisible: boolean,
): Sphere3D {
  boundaryHelper.reset();

  for (let i = 0, il = renderables.length; i < il; ++i) {
    if (onlyVisible && !renderables[i].state.visible) continue;
    if (!renderables[i].values.drawCount.ref.value) continue;

    const boundingSphere = renderables[i].values.boundingSphere.ref.value;
    if (!boundingSphere.radius) continue;

    boundaryHelper.includeSphere(boundingSphere);
  }
  boundaryHelper.finishedIncludeStep();
  for (let i = 0, il = renderables.length; i < il; ++i) {
    if (onlyVisible && !renderables[i].state.visible) continue;
    if (!renderables[i].values.drawCount.ref.value) continue;

    const boundingSphere = renderables[i].values.boundingSphere.ref.value;
    if (!boundingSphere.radius) continue;

    boundaryHelper.radiusSphere(boundingSphere);
  }

  return boundaryHelper.getSphere(boundingSphere);
}

function renderableSort(a: GraphicsRenderable, b: GraphicsRenderable) {
  const drawProgramIdA = a.getProgram('color').id;
  const drawProgramIdB = b.getProgram('color').id;
  const materialIdA = a.materialId;
  const materialIdB = b.materialId;

  if (drawProgramIdA !== drawProgramIdB) {
    // sort by program id to minimize gl state changes
    return drawProgramIdA - drawProgramIdB;
  } else if (materialIdA !== materialIdB) {
    // sort by material id to minimize gl state changes
    return materialIdA - materialIdB;
  } else {
    return a.id - b.id;
  }
}

interface Scene extends Object3D {
  readonly count: number;
  readonly renderables: ReadonlyArray<GraphicsRenderable>;
  /** all renderables as rows (`DrawTable`): the scene drawn as a group (`Renderer.renderBlended`) */
  readonly table: DrawTable;
  readonly boundingSphere: Sphere3D;
  readonly boundingSphereVisible: Sphere3D;
  readonly transparency: Transparency;

  readonly primitives: Scene.Group;
  readonly volumes: Scene.Group;

  /** Returns `true` if some visibility has changed, `false` otherwise. */
  syncVisibility: () => boolean;
  setTransparency: (transparency: Transparency) => void;
  setGlobals: (values: GlobalDefines) => void;
  update: (objects: ArrayLike<GraphicsRenderObject> | undefined, keepBoundingSphere: boolean) => void;
  /**
   * `update` of the renderables whose values changed since the last `update` of all or of the changed ones (a change of
   * a few of many objects, e.g. markers on hover), `keepBoundingSphere` as for `update`
   */
  updateChanged: (keepBoundingSphere: boolean) => void;
  add: (o: GraphicsRenderObject) => void; // GraphicsRenderable
  remove: (o: GraphicsRenderObject) => void;
  commit: (maxTimeMs?: number) => boolean;
  readonly needsCommit: boolean;
  readonly commitQueueSize: number;
  has: (o: GraphicsRenderObject) => boolean;
  /** the renderable of a committed object */
  get: (o: GraphicsRenderObject) => GraphicsRenderable | undefined;
  clear: () => void;
  forEach: (callbackFn: (value: GraphicsRenderable, key: GraphicsRenderObject) => void) => void;
  /**
   * Draws with `ctx` from now on: disposes the renderables (while the previous context is alive) and creates them
   * again on `ctx`, with textures loaded from data transferred (`transferTextures`), queued adds stay queued; those
   * with rendered or computed textures of the previous context (e.g. GPU surfaces) are dropped and returned, their
   * producers create them again
   */
  setContext: (ctx: GraphicsContext) => GraphicsRenderObject[];
  /** Marker average of primitive and direct volume renderables */
  readonly markerAverage: number;
  /** Emissive average of primitive renderables */
  readonly emissiveAverage: number;
  /** Wiggle average of primitive renderables */
  readonly wiggleAverage: number;
  /** Opacity average of primitive renderables */
  readonly opacityAverage: number;
  /** Transparency minimum, excluding fully opaque, of primitive renderables */
  readonly transparencyMin: number;
  /** Is `true` if any primitive renderable (possibly) has any opaque part */
  readonly hasOpaque: boolean;
  /** Is `true` if any primitive renderable has animation enabled */
  readonly hasAnimation: boolean;
}

namespace Scene {
  export interface Group extends Object3D {
    readonly renderables: ReadonlyArray<GraphicsRenderable>;
    /** the renderables as rows (`DrawTable`), in their order; not for groups only planned with (no scene's) */
    readonly table?: DrawTable;
  }

  export function create(
    _ctx: GraphicsContext,
    transparency: Transparency = 'blended',
    globals: GlobalDefines = { dColorMarker: true, dLightCount: 1 },
  ): Scene {
    let ctx = _ctx;
    const renderableMap = new Map<GraphicsRenderObject, GraphicsRenderable>();
    const renderables: GraphicsRenderable[] = [];
    const boundingSphere = Sphere3D();
    const boundingSphereVisible = Sphere3D();

    const primitives: GraphicsRenderable[] = [];
    const volumes: GraphicsRenderable[] = [];
    const renderablesTable = createDrawTable(renderables);
    const primitivesTable = createDrawTable(primitives);
    const volumesTable = createDrawTable(volumes);

    let boundingSphereDirty = true;
    let boundingSphereVisibleDirty = true;

    let markerAverageDirty = true;
    let emissiveAverageDirty = true;
    let wiggleAverageDirty = true;
    let opacityAverageDirty = true;
    let transparencyMinDirty = true;
    let hasOpaqueDirty = true;
    let hasAnimationDirty = true;

    let markerAverage = 0;
    let emissiveAverage = 0;
    let wiggleAverage = 0;
    let opacityAverage = 0;
    let transparencyMin = 0;
    let hasOpaque = false;
    let hasAnimation = false;

    const object3d = Object3D.create();
    const { view, position, direction, up } = object3d;

    function add(o: GraphicsRenderObject) {
      if (!renderableMap.has(o)) {
        const renderable = createRenderable(ctx, o, transparency, globals);
        renderables.push(renderable);
        renderablesTable.add(renderable);
        if (o.type === 'direct-volume') {
          volumes.push(renderable);
          volumesTable.add(renderable);
        } else {
          primitives.push(renderable);
          primitivesTable.add(renderable);
        }
        renderableMap.set(o, renderable);
        boundingSphereDirty = true;
        boundingSphereVisibleDirty = true;
      } else {
        console.warn(`RenderObject with id '${o.id}' already present`);
      }
    }

    function remove(o: GraphicsRenderObject) {
      const renderable = renderableMap.get(o);
      if (renderable) {
        renderable.dispose();
        arraySetRemove(renderables, renderable);
        arraySetRemove(primitives, renderable);
        arraySetRemove(volumes, renderable);
        renderablesTable.remove(renderable);
        primitivesTable.remove(renderable);
        volumesTable.remove(renderable);
        renderableMap.delete(o);
        boundingSphereDirty = true;
        boundingSphereVisibleDirty = true;
      }
    }

    const commitBulkSize = 100;
    function commit(maxTimeMs: number) {
      const start = now();

      let i = 0;

      while (true) {
        const o = commitQueue.tryGetRemove();
        if (!o) break;
        remove(o);
        if (++i % commitBulkSize === 0 && now() - start > maxTimeMs) return false;
      }

      while (true) {
        const o = commitQueue.tryGetAdd();
        if (!o) break;
        add(o);
        if (++i % commitBulkSize === 0 && now() - start > maxTimeMs) return false;
      }

      renderables.sort(renderableSort);
      renderablesTable.invalidateOrder();
      markerAverageDirty = true;
      emissiveAverageDirty = true;
      wiggleAverageDirty = true;
      opacityAverageDirty = true;
      transparencyMinDirty = true;
      hasOpaqueDirty = true;
      hasAnimationDirty = true;
      return true;
    }

    const commitQueue = new CommitQueue();

    function clear() {
      for (let i = 0, il = renderables.length; i < il; ++i) {
        renderablesTable.remove(renderables[i]);
        primitivesTable.remove(renderables[i]);
        volumesTable.remove(renderables[i]);
        renderables[i].dispose();
      }
      renderables.length = 0;
      primitives.length = 0;
      volumes.length = 0;
      renderableMap.clear();
      boundingSphereDirty = true;
      boundingSphereVisibleDirty = true;
    }

    function setContext(value: GraphicsContext) {
      const dropped: GraphicsRenderObject[] = [];
      const present = Array.from(renderableMap.keys());
      const queued: GraphicsRenderObject[] = [];
      while (true) {
        const o = commitQueue.tryGetAdd();
        if (!o) break;
        queued.push(o);
      }
      clear();
      ctx = value;
      // drawn again right away, as before (scenes that are not committed per frame, e.g. the helpers', stay complete);
      // queued adds stay queued, as do queued removes
      for (const o of present) {
        if (transferTextures(o, value)) add(o);
        else dropped.push(o);
      }
      for (const o of queued) {
        if (transferTextures(o, value)) commitQueue.add(o);
        else dropped.push(o);
      }
      renderables.sort(renderableSort);
      renderablesTable.invalidateOrder();
      markerAverageDirty = true;
      emissiveAverageDirty = true;
      wiggleAverageDirty = true;
      opacityAverageDirty = true;
      transparencyMinDirty = true;
      hasOpaqueDirty = true;
      hasAnimationDirty = true;
      return dropped;
    }

    /** the visibility of the renderables last synced (`DrawTable.visibleVersion` of all renderables) */
    let visibleVersion = -1;

    /** after renderables were updated: what depends on their values computed again when read */
    function updated(keepBoundingSphere: boolean) {
      if (!keepBoundingSphere) {
        boundingSphereDirty = true;
        boundingSphereVisibleDirty = true;
      } else {
        syncVisibility();
      }
      markerAverageDirty = true;
      emissiveAverageDirty = true;
      wiggleAverageDirty = true;
      opacityAverageDirty = true;
      transparencyMinDirty = true;
      hasOpaqueDirty = true;
      hasAnimationDirty = true;
    }

    function syncVisibility() {
      renderablesTable.refresh();
      if (renderablesTable.visibleVersion !== visibleVersion) {
        boundingSphereVisibleDirty = true;
        markerAverageDirty = true;
        emissiveAverageDirty = true;
        wiggleAverageDirty = true;
        opacityAverageDirty = true;
        transparencyMinDirty = true;
        hasOpaqueDirty = true;
        hasAnimationDirty = true;
        visibleVersion = renderablesTable.visibleVersion;
        return true;
      } else {
        return false;
      }
    }

    /**
     * the visible objects' slots in `t`, its rows refreshed: the averages below read its contribution columns (written
     * with the rows) instead of the objects' values
     */
    function visibleSlotsOf(t: DrawTable, slots: number[]) {
      t.refresh();
      const { order, bits, count } = t;
      slots.length = 0;
      for (let i = 0; i < count; ++i) {
        const slot = order[i];
        if (bits[slot * DrawRow.Stride + DrawRow.State] & DrawState.Visible) slots.push(slot);
      }
      return slots;
    }
    const visibleSlots: number[] = [];
    const visibleVolumeSlots: number[] = [];

    function visiblePrimitiveSlots() {
      return visibleSlotsOf(primitivesTable, visibleSlots);
    }

    function average(column: Float32Array) {
      const slots = visiblePrimitiveSlots();
      if (slots.length === 0) return 0;
      let sum = 0;
      for (let i = 0, il = slots.length; i < il; ++i) sum += column[slots[i]];
      return sum / slots.length;
    }

    /** of the primitives and the direct volumes: both are drawn by the marking pass */
    function calculateMarkerAverage() {
      const slots = visiblePrimitiveSlots();
      const volumeSlots = visibleSlotsOf(volumesTable, visibleVolumeSlots);
      const count = slots.length + volumeSlots.length;
      if (count === 0) return 0;
      let sum = 0;
      const marker = primitivesTable.contributions.marker;
      for (let i = 0, il = slots.length; i < il; ++i) sum += marker[slots[i]];
      const volumeMarker = volumesTable.contributions.marker;
      for (let i = 0, il = volumeSlots.length; i < il; ++i) sum += volumeMarker[volumeSlots[i]];
      return sum / count;
    }

    function calculateEmissiveAverage() {
      return average(primitivesTable.contributions.emissive);
    }

    function calculateWiggleAverage() {
      return average(primitivesTable.contributions.wiggle);
    }

    function calculateOpacityAverage() {
      return average(primitivesTable.contributions.opacity);
    }

    /** exclude fully opaque parts */
    function calculateTransparencyMin() {
      const slots = visiblePrimitiveSlots();
      const column = primitivesTable.contributions.transparencyMin;
      let transparencyMin = 1;
      for (let i = 0, il = slots.length; i < il; ++i) transparencyMin = Math.min(transparencyMin, column[slots[i]]);
      return transparencyMin;
    }

    function calculateHasOpaque() {
      const slots = visiblePrimitiveSlots();
      const { bits, contributions } = primitivesTable;
      for (let i = 0, il = slots.length; i < il; ++i) {
        const slot = slots[i];
        if (contributions.opaque[slot] || bits[slot * DrawRow.Stride + DrawRow.State] & DrawState.Opaque) return true;
      }
      return false;
    }

    function calculateHasAnimation() {
      const slots = visiblePrimitiveSlots();
      const { bits } = primitivesTable;
      for (let i = 0, il = slots.length; i < il; ++i) {
        if (bits[slots[i] * DrawRow.Stride + DrawRow.Facts] & DrawFact.Animated) return true;
      }
      return false;
    }

    return {
      view,
      position,
      direction,
      up,

      renderables,
      table: renderablesTable,
      primitives: { view, position, direction, up, renderables: primitives, table: primitivesTable },
      volumes: { view, position, direction, up, renderables: volumes, table: volumesTable },

      syncVisibility,
      setTransparency: (value: Transparency) => {
        transparency = value;
        for (let i = 0, il = renderables.length; i < il; ++i) {
          renderables[i].setTransparency(value);
        }
        // the programs changed
        renderablesTable.invalidate();
        primitivesTable.invalidate();
        volumesTable.invalidate();
      },
      setGlobals: (values: GlobalDefines) => {
        globals = values;
        for (const r of renderables) {
          ValueCell.updateIfChanged(r.values.dColorMarker, values.dColorMarker);
          ValueCell.updateIfChanged(r.values.dLightCount, values.dLightCount);
          r.update();
        }
      },
      update(objects, keepBoundingSphere) {
        if (objects) {
          for (let i = 0, il = objects.length; i < il; ++i) {
            const renderable = renderableMap.get(objects[i]);
            if (!renderable) continue;
            // updated: not again by the next `updateChanged` (a hover after a recolor of all), unless changed meanwhile
            const slot = renderablesTable.getSlot(renderable);
            if (slot >= 0) renderablesTable.dropChanged(slot);
            renderable.update();
          }
        } else {
          // all updated: so are the changed ones
          renderablesTable.takeChanged();
          for (let i = 0, il = renderables.length; i < il; ++i) {
            renderables[i].update();
          }
        }
        updated(keepBoundingSphere);
      },
      updateChanged(keepBoundingSphere) {
        const slots = renderablesTable.takeChanged();
        const { objects } = renderablesTable;
        for (let i = 0, il = slots.length; i < il; ++i) objects[slots[i]]?.update();
        updated(keepBoundingSphere);
      },
      add: (o: GraphicsRenderObject) => commitQueue.add(o),
      remove: (o: GraphicsRenderObject) => commitQueue.remove(o),
      commit: (maxTime = Number.MAX_VALUE) => commit(maxTime),
      get commitQueueSize() {
        return commitQueue.size;
      },
      get needsCommit() {
        return !commitQueue.isEmpty;
      },
      has: (o: GraphicsRenderObject) => {
        return renderableMap.has(o);
      },
      get: (o: GraphicsRenderObject) => renderableMap.get(o),
      clear,
      forEach: (callbackFn: (value: GraphicsRenderable, key: GraphicsRenderObject) => void) => {
        renderableMap.forEach(callbackFn);
      },
      setContext,
      get count() {
        return renderables.length;
      },
      get boundingSphere() {
        if (boundingSphereDirty) {
          calculateBoundingSphere(renderables, boundingSphere, false);
          boundingSphereDirty = false;
        }
        return boundingSphere;
      },
      get boundingSphereVisible() {
        if (boundingSphereVisibleDirty) {
          calculateBoundingSphere(renderables, boundingSphereVisible, true);
          boundingSphereVisibleDirty = false;
        }
        return boundingSphereVisible;
      },
      get transparency() {
        return transparency;
      },
      get markerAverage() {
        if (markerAverageDirty) {
          markerAverage = calculateMarkerAverage();
          markerAverageDirty = false;
        }
        return markerAverage;
      },
      get emissiveAverage() {
        if (emissiveAverageDirty) {
          emissiveAverage = calculateEmissiveAverage();
          emissiveAverageDirty = false;
        }
        return emissiveAverage;
      },
      get wiggleAverage() {
        if (wiggleAverageDirty) {
          wiggleAverage = calculateWiggleAverage();
          wiggleAverageDirty = false;
        }
        return wiggleAverage;
      },
      get opacityAverage() {
        if (opacityAverageDirty) {
          opacityAverage = calculateOpacityAverage();
          opacityAverageDirty = false;
        }
        return opacityAverage;
      },
      get transparencyMin() {
        if (transparencyMinDirty) {
          transparencyMin = calculateTransparencyMin();
          transparencyMinDirty = false;
        }
        return transparencyMin;
      },
      get hasOpaque() {
        if (hasOpaqueDirty) {
          hasOpaque = calculateHasOpaque();
          hasOpaqueDirty = false;
        }
        return hasOpaque;
      },
      get hasAnimation() {
        if (hasAnimationDirty) {
          hasAnimation = calculateHasAnimation();
          hasAnimationDirty = false;
        }
        return hasAnimation;
      },
    };
  }
}

export { Scene };
