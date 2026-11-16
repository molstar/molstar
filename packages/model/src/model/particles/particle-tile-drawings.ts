/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Clip } from '@molstar/core/util/clip';
import type { ParticleList } from './particle-list.js';
import { ParticleTileView } from './particle-tile-view.js';
import type { ParticleTileDrawing } from './particle-tile-visibility.js';
import type { ParticleTileSet } from './particle-tiles.js';

// What draws a tile set's views (`ParticleTileDrawing`: a representation each, its clip objects and whether opaque):
// what a view by the camera picks its tiles by (`ParticleTileVisibilityOptions.drawings`). Set by the representations
// as they update, a clipping gizmo's drag too (an update of the representation alone, not the state's), so the picks
// follow a cut as it moves; told per set (`subscribe`) to who picks again.

/** what a representation draws with: its clip objects, its alpha and x-ray (a representation's props) */
export interface ParticleTileDrawingProps {
  readonly clip?: Clip.Props;
  readonly alpha?: number;
  readonly xrayShaded?: unknown;
}

/**
 * a set's drawings, a version bumped as they change (a pick's key: picked again), when they last changed
 * (`performance.now()`'s: a cut still moving while within `ParticleTileDrawings.SettleMs`)
 */
export interface ParticleTileDrawingsOf {
  readonly list: readonly ParticleTileDrawing[];
  readonly version: number;
  readonly changedAt: number;
}

interface Entry {
  readonly drawing: ParticleTileDrawing;
  readonly key: string;
}

interface Drawings {
  readonly byOwner: Map<object, Entry>;
  list: ParticleTileDrawing[];
  version: number;
  changedAt: number;
  readonly listeners: Set<() => void>;
}

const bySet = new WeakMap<ParticleTileSet, Drawings>();
/** per owner the set it draws, to let go of where it draws another or none */
const drawnSets = new WeakMap<object, ParticleTileSet>();

function getDrawings(set: ParticleTileSet): Drawings {
  let d = bySet.get(set);
  if (!d)
    bySet.set(set, (d = { byOwner: new Map(), list: [], version: 0, changedAt: -Infinity, listeners: new Set() }));
  return d;
}

function changed(d: Drawings) {
  d.list = [...d.byOwner.values()].map((e) => e.drawing);
  d.version += 1;
  d.changedAt = performance.now();
  for (const l of d.listeners) l();
}

export namespace ParticleTileDrawings {
  /**
   * a set's drawings changed within this many ms: still changing (a cut dragged), its picks the quicker ones (a coarser
   * occlusion), picked again as they settle
   */
  export const SettleMs = 250;

  /** whether `set`'s drawings changed within `SettleMs` of `now` (a cut being moved) */
  export function isChanging(set: ParticleTileSet, now = performance.now()) {
    return now - getDrawings(set).changedAt < SettleMs;
  }

  /**
   * How `owner` (a representation) draws `list` (a tile view's; else, or `props` none: it draws none, let go of) with
   * `props`: its clip objects (any shapes) and whether opaque (alpha 1, no x-ray). Its set's listeners told where that
   * changed (a moved clip object), not where the same again
   */
  export function set(owner: object, list: ParticleList | undefined, props: ParticleTileDrawingProps | undefined) {
    const set = list ? ParticleTileView.of(list)?.set : undefined;
    const last = drawnSets.get(owner);
    if (last && (last !== set || !props)) {
      const d = getDrawings(last);
      if (d.byOwner.delete(owner)) changed(d);
      drawnSets.delete(owner);
    }
    if (!set || !props) return;
    drawnSets.set(owner, set);
    const objects = props.clip?.objects ?? [];
    const opaque = (props.alpha ?? 1) >= 1 && !props.xrayShaded;
    const key = `${opaque},${props.clip?.variant},${JSON.stringify(objects)}`;
    const d = getDrawings(set);
    if (d.byOwner.get(owner)?.key === key) return;
    const clip = objects.length > 0 ? Clip.getClip(props.clip!).objects : undefined;
    d.byOwner.set(owner, { drawing: { clip, opaque }, key });
    changed(d);
  }

  /** `set`'s drawings now (none: drawn whole and opaque) */
  export function get(set: ParticleTileSet): ParticleTileDrawingsOf {
    return getDrawings(set);
  }

  /** `listener` told as `set`'s drawings change; the returned function unsubscribes */
  export function subscribe(set: ParticleTileSet, listener: () => void): () => void {
    const d = getDrawings(set);
    d.listeners.add(listener);
    return () => d.listeners.delete(listener);
  }
}
