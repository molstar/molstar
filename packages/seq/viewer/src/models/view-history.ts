/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { hash2, hashString } from '@molstar/core/data/util/hash-functions';
import { type ReadonlyState, State } from '@molstar/core/reactive';
import type { Location } from '@molstar/seqx-core/query';
import type { BandViewport } from '@molstar/seqx-core/render';
import {
  placeOfUnit,
  type ResolvedBlock,
  type ResolvedSection,
  type SystemPosition,
  unitOfPlace,
} from '@molstar/seqx-core/view';
import { axisKey, axisTotal, positionAt } from './axis.js';
import type { ViewBands } from './view-bands.js';
import { centerOffset, DetachedPrefix, type Group, type SectionState } from './view-groups.js';

// The view's history of jumps (_spec/README.md "ViewModel", History): a step is every group's range with its axis's
// signature, every band's top row by key, every grid's own y window; restored onto the view as it is now (by position
// onto a changed axis, a follower detached or attached again as it was); the steps back and forward, a step to where
// the view is now skipped.

/**
 * A group's range in a step of the history (_spec/README.md "ViewModel", History): its offsets with its axis's
 * signature (system, pieces, direction), its centre as a position and its members, so a step over a changed axis is
 * restored by position (the centre and the zoom), never dropped.
 */
interface Kept {
  readonly start: number;
  readonly end: number;
  readonly fitted: boolean;
  readonly sig?: string;
  readonly center?: number;
  readonly members?: readonly string[];
  readonly scale?: number;
  /**
   * On a section of blocks, its ends as `{ system, position }` (never units): a step over a changed block list (an
   * `include`, an order, a scale) is restored by them, an end on a block no longer shown at the nearest shown.
   */
  readonly ends?: readonly [SystemPosition, SystemPosition];
}

/**
 * A band as a step of the history keeps it (_spec/README.md "ViewModel", History): its top row by key (a row's key, and
 * the fraction of it scrolled past) with its row zoom and folds, so a step back lands on the same rows whatever was
 * folded, hidden or reordered since; `bandTop` as it was, where no layout names the rows (`BandRows`).
 */
export interface KeptBand {
  readonly band: BandViewport;
  readonly key?: string;
  readonly offset?: number;
}

/**
 * The leading viewports at one time: each group's range by key (a detached follower's among them), the bands, and the
 * grids' own y windows (by section; a grid absent shows all of its y).
 */
export interface Snapshot {
  readonly groups: ReadonlyMap<string, Kept>;
  readonly bands: ReadonlyMap<string, KeptBand>;
  readonly ys: ReadonlyMap<string, { readonly start: number; readonly end: number }>;
  /** The Location the jump that left the view here named, while no hand moved it since (`ViewModel.lastJump`). */
  readonly location?: Location;
}

/** What a jump was (_spec/README.md "ViewModel", History, The last jump). */
export type JumpKind = 'focus' | 'reveal' | 'back' | 'forward' | 'fit' | 'zoom' | 'move';

/**
 * The last jump (`ViewModel.lastJump`): its count (`jumps`), its kind, the Location it named (a focus's or a reveal's;
 * a step's when the step kept one) and the duration the viewer animates it over.
 */
export interface LastJump {
  readonly n: number;
  readonly kind: JumpKind;
  readonly location?: Location;
  readonly durationMs: number;
}

/** A step of the view's history as `snapshot` gives it: opaque, only `record` reads it. */
export interface ViewSnapshot {
  readonly __viewSnapshot?: never;
}

/** What a jump did: moved the view, found it already there, or could not move. */
export type JumpOutcome = 'moved' | 'same' | 'none';

/** The view's history of jumps: how many steps back and forward there are. */
export interface ViewHistory {
  readonly back: number;
  readonly forward: number;
}

const sameBand = (a: KeptBand | undefined, b: KeptBand | undefined) =>
  !!a &&
  !!b &&
  a.band.rowScale === b.band.rowScale &&
  a.band.folds === b.band.folds &&
  (a.key !== undefined && b.key !== undefined
    ? a.key === b.key && a.offset === b.offset
    : a.band.bandTop === b.band.bandTop);

/** Two steps are the same view: the same groups at the same offsets, compared only on equal axes. */
export function snapshotsEqual(a: Snapshot, b: Snapshot): boolean {
  if (a.groups.size !== b.groups.size || a.bands.size !== b.bands.size) return false;
  for (const [k, v] of a.groups) {
    const w = b.groups.get(k);
    if (!w || w.sig !== v.sig || w.start !== v.start || w.end !== v.end || w.fitted !== v.fitted || w.scale !== v.scale)
      return false;
  }
  for (const [k, v] of a.bands) if (!sameBand(v, b.bands.get(k))) return false;
  if (a.ys.size !== b.ys.size) return false;
  for (const [k, v] of a.ys) {
    const w = b.ys.get(k);
    if (!w || w.start !== v.start || w.end !== v.end) return false;
  }
  return true;
}

/** Each group's range now, with its axis's signature, its centre as a position and its members. */
export function snapshotOf(
  groups: ReadonlyMap<string, Group>,
  sections: ReadonlyMap<string, SectionState>,
  bands: ViewBands,
): Snapshot {
  const out = new Map<string, Kept>();
  for (const g of groups.values()) {
    const center = g.axis.ranges ? positionAt(g.axis, (g.start + g.end) / 2) : undefined;
    const blocks = blocksOf(g, sections);
    const a = blocks && placeOfUnit(blocks, g.start, 'start');
    const b = blocks && placeOfUnit(blocks, g.end, 'end');
    out.set(g.key, {
      start: g.start,
      end: g.end,
      fitted: g.fitted,
      sig: sigOf(g, sections),
      ...(center !== undefined ? { center } : {}),
      ...(a && b ? { ends: [a, b] as const } : {}),
      members: g.members,
      ...(g.scale !== undefined ? { scale: g.scale } : {}),
    });
  }
  // the bands with the windows, by their top rows' keys: back after a fold or a scroll lands on the same rows
  const kept = new Map<string, KeptBand>();
  for (const [track, band] of bands.value) {
    const at = bands.rows?.keyAt(track, band);
    kept.set(track, at ? { band, key: at.key, offset: at.offset } : { band });
  }
  // a grid's own y window (`follow_y: none`); a y that follows a section is kept by that section's group
  const ys = new Map<string, { start: number; end: number }>();
  for (const [address, s] of sections) if (s.y) ys.set(address, { ...s.y });
  return { groups: out, bands: kept, ys };
}

/**
 * A group's kept range in a step: by its key, else (its axis changed since: a reverse, a collapse) by a member; on an
 * axis of another signature, restored by position: the centre and the zoom, the centre clamped to what is drawn (a
 * centre now hidden goes to its break). Undefined when the step has none of it.
 */
function keptOf(snap: Snapshot, g: Group, sections: ReadonlyMap<string, SectionState>): Kept | undefined {
  let v = snap.groups.get(g.key);
  if (!v)
    for (const w of snap.groups.values())
      if (w.members?.some((m) => g.members.includes(m))) {
        v = w;
        break;
      }
  const sig = sigOf(g, sections);
  if (!v || v.sig === undefined || v.sig === sig) return v;
  // the letter size is kept within one layout only (a text's, across a reverse or a collapse)
  const scale = v.scale !== undefined && layoutOfSig(v.sig) === g.layout ? { scale: v.scale } : {};
  const total = axisTotal(g.axis);
  if (v.fitted || total === undefined) return { start: 0, end: total ?? v.end, fitted: true, sig, ...scale };
  // a section of blocks: its ends by system and position onto the blocks now (clamped to what is drawn)
  const blocks = v.ends && blocksOf(g, sections);
  if (blocks && v.ends) {
    const a = unitOfPlace(blocks, v.ends[0]);
    const b = unitOfPlace(blocks, v.ends[1], a);
    // a reversed block since maps the ends mirrored (b before a): the range is their hull either way
    if (a !== undefined && b !== undefined && a !== b)
      return { start: Math.min(a, b), end: Math.max(a, b), fitted: false, sig, ...scale };
  }
  const c = v.center === undefined ? undefined : centerOffset(g.axis, v.center);
  if (c === undefined) return { start: 0, end: total, fitted: true, sig, ...scale };
  // the centre and the zoom: a window's middle, a text's lines in view (its anchor and its estimated end) alike
  const half = (v.end - v.start) / 2;
  return { start: c - half, end: c + half, fitted: false, sig, ...scale };
}

/**
 * What restoring a snapshot would give: the groups there now, at its ranges (by position onto a changed axis); a
 * follower detached since, attached; one detached then and attached since, detached again (a section gone since is
 * left out). `ownY`: whether a section is a grid with a y window of its own.
 */
export function restoredOf(
  snap: Snapshot,
  groups: ReadonlyMap<string, Group>,
  sections: ReadonlyMap<string, SectionState>,
  ownY: (section: string) => boolean,
): Snapshot {
  const out = new Map<string, Kept>();
  for (const g of groups.values()) {
    const v = keptOf(snap, g, sections);
    if (v) out.set(g.key, { ...v, sig: sigOf(g, sections) });
    else if (!g.key.startsWith(DetachedPrefix))
      out.set(g.key, {
        start: g.start,
        end: g.end,
        fitted: g.fitted,
        sig: sigOf(g, sections),
        ...(g.scale !== undefined ? { scale: g.scale } : {}),
      });
  }
  for (const [key, v] of snap.groups)
    if (key.startsWith(DetachedPrefix) && !out.has(key) && detachable(sections, key)) out.set(key, v);
  const ys = new Map<string, { start: number; end: number }>();
  for (const [address, y] of snap.ys) if (ownY(address)) ys.set(address, y);
  return { groups: out, bands: snap.bands, ys };
}

/**
 * The groups take the snapshot's ranges (the caller publishes them in one go); a follower detached since it was taken
 * is attached again, one detached then is detached again; each band at its kept top row's place; each grid's own y as
 * kept.
 */
export function restoreInto(
  snap: Snapshot,
  groups: Map<string, Group>,
  sections: ReadonlyMap<string, SectionState>,
  bands: ViewBands,
  ownY: (section: string) => boolean,
): void {
  // each band at its kept top row's place under its kept folds and zoom, clamped; a band it had none of has none
  const map = new Map<string, BandViewport>();
  for (const [track, k] of snap.bands) {
    const top = k.key !== undefined ? bands.rows?.topOf(track, k.band, k.key, k.offset ?? 0) : undefined;
    map.set(track, top === undefined ? k.band : { ...k.band, bandTop: top });
  }
  if (map.size || bands.value.size) bands.replace(map);
  for (const [key, g] of [...groups]) {
    const v = keptOf(snap, g, sections);
    if (v) {
      g.start = v.start;
      g.end = v.end;
      g.fitted = v.fitted;
      if (v.scale !== undefined) g.scale = v.scale;
      else delete g.scale;
    } else if (key.startsWith(DetachedPrefix)) {
      groups.delete(key);
      for (const s of sections.values()) if (s.group === key) s.group = undefined;
    }
  }
  // a follower detached when the snapshot was taken is detached again, at its range then
  for (const [key, v] of snap.groups) {
    const s = !groups.has(key) && key.startsWith(DetachedPrefix) ? detachable(sections, key) : undefined;
    if (!s) continue;
    groups.set(key, {
      key,
      axis: s.section.axis,
      layout: s.section.layout,
      start: v.start,
      end: v.end,
      fitted: v.fitted,
      ...(v.scale !== undefined ? { scale: v.scale } : {}),
      placedWidth: s.width,
      members: [s.section.address],
    });
    s.group = key;
  }
  // each grid's own y as kept (none kept: all of its y)
  for (const [address, s] of sections) {
    if (!ownY(address)) continue;
    const y = snap.ys.get(address);
    if (y) s.y = { ...y };
    else delete s.y;
  }
}

/** The follower a detached group's key names, when it is in the view and follows. */
function detachable(sections: ReadonlyMap<string, SectionState>, key: string): SectionState | undefined {
  const s = sections.get(key.slice(DetachedPrefix.length));
  return s?.section.follow ? s : undefined;
}

/**
 * A group's signature in a step of the history: its axis's and its layout, and on a section of blocks its blocks'
 * systems, directions, scales and offsets (a step over any is restored by position), hashed: a few bytes a step
 * whatever the blocks (3,000 contigs were ~150 kB a step written out), memoized by the blocks list.
 */
function sigOf(g: Group, sections: ReadonlyMap<string, SectionState>): string {
  const blocks = blocksOf(g, sections)?.blocks;
  return `${axisKey(g.axis)}#${g.layout}${blocks ? `#${blocksSig(blocks)}` : ''}`;
}

const BlocksSigs = new WeakMap<readonly ResolvedBlock[], string>();

/** A blocks list's signature: its count and two 32-bit hashes of each block's system, direction, scale and offset. */
export function blocksSig(blocks: readonly ResolvedBlock[]): string {
  const known = BlocksSigs.get(blocks);
  if (known !== undefined) return known;
  let h1 = 0;
  let h2 = 0x9e3779b9 | 0;
  for (const b of blocks) {
    const s = hashString(b.axis.system);
    const d = b.axis.dir === -1 ? 1 : 0;
    const o = Math.round((b.offset ?? 0) * 1024) | 0;
    const k = Math.round((b.scale ?? 1) * 1e6) | 0;
    h1 = hash2(hash2(h1, s), hash2(hash2(d, k), o)) | 0;
    h2 = (Math.imul(h2 ^ s, 0x01000193) + Math.imul(o ^ (k << 1) ^ d, 0x5bd1e995)) | 0;
  }
  const sig = `|${blocks.length}:${(h1 >>> 0).toString(36)}:${(h2 >>> 0).toString(36)}`;
  BlocksSigs.set(blocks, sig);
  return sig;
}

/** A group's section of blocks (its first member's), when it is one. */
function blocksOf(g: Group, sections: ReadonlyMap<string, SectionState>): ResolvedSection | undefined {
  const s = sections.get(g.members[0])?.section;
  return s?.blocks?.length ? s : undefined;
}

function layoutOfSig(sig: string): string {
  // the layout after the axis's key, before a blocks' key
  const parts = sig.split('#');
  return parts.length > 2 && parts[parts.length - 1].startsWith('|')
    ? parts[parts.length - 2]
    : parts[parts.length - 1];
}

/** The steps back and forward: what `jump` and `record` push, `back` and `forward` take. */
export class JumpSteps {
  private readonly _back: Snapshot[] = [];
  private readonly _forward: Snapshot[] = [];
  private readonly _history = new State<ViewHistory>(
    { back: 0, forward: 0 },
    (a, b) => a.back === b.back && a.forward === b.forward,
  );
  /** How many steps back and forward there are. */
  readonly history: ReadonlyState<ViewHistory> = this._history.asReadonly;

  get canBack(): boolean {
    return this._back.length > 0;
  }

  get canForward(): boolean {
    return this._forward.length > 0;
  }

  /** Pushes a step back, at most `limit` kept, the forward steps dropped. */
  push(kept: Snapshot, limit: number): void {
    this._back.push(kept);
    const max = Math.max(0, limit);
    if (this._back.length > max) this._back.splice(0, this._back.length - max);
    this._forward.length = 0;
    this._tell();
  }

  /** Forgets every step. */
  clear(): void {
    this._back.length = 0;
    this._forward.length = 0;
    this._tell();
  }

  /**
   * A step back (or forward): the first step whose restoring is not the view `now` (one to where the view is now, a
   * group gone, a hand that returned, is skipped) is restored, `now` pushed the other way. False when there is none.
   */
  step(
    forward: boolean,
    now: Snapshot,
    restored: (snap: Snapshot) => Snapshot,
    restore: (snap: Snapshot) => void,
  ): boolean {
    const from = forward ? this._forward : this._back;
    const to = forward ? this._back : this._forward;
    let target: Snapshot | undefined;
    while (from.length) {
      const t = from.pop()!;
      if (!snapshotsEqual(restored(t), now)) {
        target = t;
        break;
      }
    }
    if (!target) {
      this._tell();
      return false;
    }
    to.push(now);
    restore(target);
    this._tell();
    return true;
  }

  private _tell(): void {
    this._history.set({ back: this._back.length, forward: this._forward.length });
  }
}
