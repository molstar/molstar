/**
 * Copyright (c) 2018-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

// A canvas's marking: loci marked (highlight, select) buffered until the next render, then applied to the
// representations they may change (`MarkIndex`) and the camera axes; and the loci of a picking id

import type { PickingId } from '@molstar/graphics/geo/geometry/picking';
import { Representation } from '@molstar/graphics/repr/representation';
import type { MarkerAction } from '@molstar/core/util/marker-action';
import { EmptyLoci, isEmptyLoci, type Loci } from '@molstar/model/model/loci';
import type { Canvas3DParts } from './parts.js';
import type { Canvas3DRepresentations } from './representations.js';

/** What the marking reads of the canvas, and calls back */
export interface Canvas3DMarksHost {
  readonly parts: Canvas3DParts;
  readonly reprs: Canvas3DRepresentations;
  /** the canvas has work (wakes its loop) */
  busy: () => void;
  /** marks were applied (`Canvas3D.interaction`'s event) */
  marked: () => void;
}

export interface Canvas3DMarks {
  /** The loci of `pickingId`, of the camera axes or of the representation that has it (`Canvas3D.getLoci`) */
  getLoci: (pickingId: PickingId | undefined) => Representation.Loci;
  /** Buffers a mark for the next render (`Canvas3D.mark`) */
  mark: (reprLoci: Representation.Loci, action: MarkerAction) => void;
  /** Applies the buffered marks: whether any changed what is drawn */
  resolve: () => boolean;
  /** marks buffered, not yet applied */
  readonly pending: boolean;
  /**
   * whether the last `resolve` changed the marking of the helpers or the overlay: they are drawn into the scene color, so
   * their marking can't be redrawn on its own
   */
  readonly helperUpdated: boolean;
  /** the buffered marks dropped */
  clear: () => void;
}

export function createCanvas3DMarks(host: Canvas3DMarksHost): Canvas3DMarks {
  const { parts, reprs } = host;

  function getLoci(pickingId: PickingId | undefined) {
    let loci: Loci = EmptyLoci;
    let repr: Representation.Any = Representation.Empty;
    if (pickingId) {
      const cameraHelperLoci = parts.helper.axes.getLoci(pickingId);
      if (cameraHelperLoci !== EmptyLoci) return { loci: cameraHelperLoci, repr };

      for (const _repr of reprs.keys()) {
        const _loci = _repr.getLoci(pickingId);
        if (!isEmptyLoci(_loci)) {
          if (!isEmptyLoci(loci)) {
            console.warn('found another loci, this should not happen');
          }
          loci = _loci;
          repr = _repr;
        }
      }
    }
    return { loci, repr };
  }

  const markBuffer: [reprLoci: Representation.Loci, action: MarkerAction][] = [];

  function mark(reprLoci: Representation.Loci, action: MarkerAction) {
    // NOTE: might try to optimize a case with opposite actions for the
    //       same loci. Tho this might end up being more expensive (and error prone)
    //       then just applying everything "naively".
    markBuffer.push([reprLoci, action]);
    host.busy();
  }

  let helperMarkingUpdated = false;

  function resolve() {
    let changed = false;
    helperMarkingUpdated = false;
    // every render: no array made where nothing was marked
    for (let i = 0; i < markBuffer.length; ++i) {
      changed = applyMark(markBuffer[i][0], markBuffer[i][1]) || changed;
    }
    markBuffer.length = 0;
    if (changed) {
      const { scene, helper, shaderManager } = parts;
      // the marked objects' values changed: those updated, not every object of the scene
      scene.updateChanged(true);
      helper.overlay.scene.updateChanged(true);
      for (const s of helper.screen.scenes) s.updateChanged(true);

      shaderManager.updateRequired(parts.props);
      shaderManager.finalizeRequired(true);

      host.marked();
    }
    return changed;
  }

  function applyMark(reprLoci: Representation.Loci, action: MarkerAction) {
    const { repr, loci } = reprLoci;
    const { scene } = parts;
    let changed = false;
    // the overlay is drawn into the scene color, as the helpers are
    let overlayChanged = false;
    const markRepr = (r: Representation.Any) => {
      if (!r.mark(loci, action)) return;
      reprs.markIndex.note(r, action);
      changed = true;
      if (reprs.sceneOf(r) !== scene) overlayChanged = true;
    };
    if (repr) {
      markRepr(repr);
    } else {
      reprs.markIndex.forEach(loci, action, markRepr);
    }
    const cameraChanged = parts.helper.axes.mark(loci, action);
    const helperChanged = cameraChanged || overlayChanged;
    helperMarkingUpdated = helperMarkingUpdated || helperChanged;
    return changed || helperChanged;
  }

  return {
    getLoci,
    mark,
    resolve,
    get pending() {
      return markBuffer.length > 0;
    },
    get helperUpdated() {
      return helperMarkingUpdated;
    },
    clear: () => {
      markBuffer.length = 0;
    },
  };
}
