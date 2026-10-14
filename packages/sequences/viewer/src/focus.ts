/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Focus (seqspec:interaction#focus), its view part: a section's resolved `focus` is its initial value; a patch or a
// snapshot that changes it sets it again, as a focus change the `focusViewport` behavior moves to. A focus from
// extents is the hull of their projection (no closure) on the first block of each section that it lights. Snapshots and
// their player are VC2's. Loci and `Coordinates` go through the edge-rate `Projector` to extents (IN1).

import type { ResolvedView } from '@molstar/seq-core/format';
import { Extents, type SeqSet } from '@molstar/seq-core/model';
import type { SeqLociAny } from './lineage.js';
import { BlockHulls, type Coordinates } from './projection.js';
import { event, type Origin, type SeqEvent, type SeqEventEmitter, type SeqState, state } from './state.js';
import type { Viewports } from './viewport.js';

/** A block's range [start, stop), block coordinate system indices; the block's index in its section. */
export interface FocusEntry {
  /** The first block it lights, in axis order, and its range there. */
  readonly block: number;
  readonly start: number;
  readonly stop: number;
  /** When it lights several blocks: each one's (block, start, stop), in axis order, the first included; else null. */
  readonly blocks?: Int32Array | null;
  /** What it was asked from (the caller's object); null: a resolved focus or `set`. */
  readonly loci: SeqLociAny | Coordinates | Extents | null;
}
/** Reused. `transitionMs` -1: not given (the behavior's `Viewport.AnimationMs`); 0: at once. */
export interface FocusChange {
  section: string;
  entry: FocusEntry | null;
  transitionMs: number;
  origin: Origin;
}
export interface FocusManager {
  /** By section ref. */
  readonly current: SeqState<ReadonlyMap<string, FocusEntry>>;
  readonly changed: SeqEvent<FocusChange>;
  set(section: string, entry: FocusEntry | null, opts?: { transitionMs?: number; origin?: Origin }): void;
  /** Every section that shows them (or `section`), as the hull of their projection on every block it lights; false: none does. */
  fromLoci(
    target: SeqLociAny | Coordinates,
    opts?: { section?: string; transitionMs?: number; origin?: Origin },
  ): boolean;
  /** The core form `fromLoci` converts to. */
  fromExtents(ext: Extents, opts?: { section?: string; transitionMs?: number; origin?: Origin }): boolean;
}

const hulls = BlockHulls.create();

/** A block by what it shows, not its index: its unit's or alignment's key, its piece and its direction. */
function blockKey(view: ResolvedView, set: SeqSet, b: number): string {
  const { unit, alignment, start, stop, dir } = view.blocks;
  const a = alignment[b] >= 0 ? set.alignments[alignment[b]] : null;
  const of = unit[b] >= 0 ? `u${set.strings.get(set.units.key[unit[b]])}` : `a${a ? set.strings.get(a.key) : ''}`;
  return `${of}:${start[b]}:${stop[b]}:${dir[b]}`;
}

export class Focus implements FocusManager {
  readonly current = state<ReadonlyMap<string, FocusEntry>>(new Map());
  readonly changed: SeqEventEmitter<FocusChange>;
  /** Set by the snapshot manager while it applies one: the transition its focus changes carry. */
  transition: { ms: number; origin: Origin } | null = null;
  private readonly resolved = new Map<string, string>();
  /** The entry the document's focus last gave each section (what a quiet move of its block replaces). */
  private readonly given = new Map<string, FocusEntry | null>();
  private readonly record: FocusChange = { section: '', entry: null, transitionMs: -1, origin: 'api' };
  private readonly viewports: Viewports;
  private readonly ext = Extents.create(4);

  constructor(viewports: Viewports, onError: (e: unknown) => void) {
    this.viewports = viewports;
    this.changed = event(onError);
  }

  /**
   * The current view's resolved focus per section: the first is the initial value, a later change sets it again. A
   * focus is the same while its block (what it shows: its unit or alignment, piece and direction) and range are: a
   * reassembly that moves that block to another index (a block inserted before it) changes no focus, and the entry the
   * document gave follows its block quietly.
   */
  update(view: ResolvedView | null, set: SeqSet | null) {
    const s = view?.sections;
    for (let i = 0; view && set && s && i < s.count; i++) {
      const ref = this.viewports.sectionRef(i) as string;
      const f = s.focus[i];
      const entry = f && { block: f.block - s.blockStart[i], start: f.start, stop: f.stop, loci: null };
      const key = f ? `${blockKey(view, set, f.block)}:${f.start}:${f.stop}` : '';
      const prev = this.resolved.get(ref);
      this.resolved.set(ref, key);
      const given = this.given.get(ref);
      if (prev === key) {
        // The same focus: the entry the document gave, while it stands, follows its block to its index.
        if (entry && given && given.block !== entry.block && this.current.value.get(ref) === given) {
          this.given.set(ref, entry);
          this.put(ref, entry, 'graph');
        }
        continue;
      }
      this.given.set(ref, entry);
      const t = this.transition;
      if (prev === undefined) this.put(ref, entry, 'graph');
      else this.set(ref, entry, { transitionMs: t?.ms ?? -1, origin: t?.origin ?? 'graph' });
    }
  }

  set(section: string, entry: FocusEntry | null, opts: { transitionMs?: number; origin?: Origin } = {}) {
    const origin = opts.origin ?? 'api';
    this.put(section, entry, origin);
    Object.assign(this.record, { section, entry, transitionMs: opts.transitionMs ?? -1, origin });
    this.changed.emit(this.record);
  }

  fromLoci(target: SeqLociAny | Coordinates, opts: { section?: string; transitionMs?: number; origin?: Origin } = {}) {
    this.ext.count = 0;
    this.viewports.projector.fromLoci(target, this.ext);
    return this.fromExtents(this.ext, opts, target);
  }

  fromExtents(
    ext: Extents,
    opts: { section?: string; transitionMs?: number; origin?: Origin } = {},
    loci: SeqLociAny | Coordinates | Extents = ext,
  ): boolean {
    let any = false;
    for (const ref of opts.section ? [opts.section] : this.viewports.current) {
      const n = this.viewports.hull(this.viewports.sectionIndex(ref), ext, hulls);
      if (n === 0) continue;
      const d = hulls.data;
      this.set(ref, { block: d[0], start: d[1], stop: d[2], blocks: n > 1 ? d.slice(0, 3 * n) : null, loci }, opts);
      any = true;
    }
    return any;
  }

  private put(section: string, entry: FocusEntry | null, origin: Origin) {
    const next = new Map(this.current.value);
    if (entry) next.set(section, entry);
    else next.delete(section);
    this.current.set(next, origin);
  }

  dispose() {
    this.current.complete();
    this.changed.complete();
  }
}
