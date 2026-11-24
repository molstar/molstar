/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { resolveAddress } from '@molstar/core/repr/tree';
import { PlacePoint } from '@molstar/seqx-core/render';
import type { ModifierState } from '@molstar/seqx-viewer/bindings';
import type { HitTarget } from '@molstar/seqx-viewer/models';
import type { SeqViewer } from './viewer.js';

// The zoom of a ring nobody follows is a cut (_spec/README.md "Input", The ring's zoom is a cut): a ring shows its whole axis, so its window is idle; the zoom gesture opens the chromosome under the
// pointer (`section.include` of its block: in) or shows every block again (out), one cut a wheel gesture (the deltas
// past `navigate.cutNotch` px, then nothing until the wheel rests `navigate.cutRestMs`). Over the hole the ring says
// once a gesture what it shows. A ring a section follows keeps the zoom of its window. A cut is the command the stem
// menu's "Show only" runs, one patch and one step of the view's history: nothing is drawn or kept per frame.

const Hole = 'The ring shows the whole genome: open a chromosome from its label’s menu or the square';

export class RingCut {
  /** When the last wheel of the gesture came (ms); the gesture is over after a rest. */
  private _last = Number.NEGATIVE_INFINITY;
  /** The gesture's deltas summed (negative: in). */
  private _sum = 0;
  /** The gesture has cut or said: nothing more until it rests. */
  private _done = false;
  /** The point under the pointer, reused. */
  private readonly _point = new PlacePoint();
  private readonly _at = new Float64Array(2);

  constructor(private readonly viewer: SeqViewer) {}

  /**
   * A zoom gesture over a section at a page point; `dy` the wheel's vertical delta in px (negative: in). True when the
   * section is an idle ring, whose zoom this is (the event is taken whatever it did); false: the zoom is the window's.
   */
  zoom(
    section: string,
    mods: ModifierState,
    target: HitTarget,
    dy: number,
    x: number,
    y: number,
    pinch = false,
  ): boolean {
    const v = this.viewer;
    if (v.bindings.resolve(pinch ? 'pinch' : 'wheel', mods, target) !== 'zoom') return false;
    const s = v.view.resolved.value?.sections.find((x) => x.address === section);
    if (!s?.blocks || !v.placementOf(section)?.warpOf || v.view.isFollowing(section)) return false;
    // a ring a section follows moves its window (its followers' windows come from it)
    if (v.view.resolved.value?.sections.some((x) => x.follow?.leader === section)) return false;
    const { cutNotch, cutRestMs } = v.config.navigate;
    const now = performance.now();
    if (now - this._last > cutRestMs) {
      this._sum = 0;
      this._done = false;
    }
    this._last = now;
    if (this._done) return true;
    this._sum += dy;
    if (Math.abs(this._sum) < cutNotch) return true;
    const include = includeOf(v, section);
    // the block whose band is under the pointer; none over the hole, the labels or a gap
    const p = v.placementOf(section)!;
    this._at[0] = x;
    this._at[1] = y;
    const block = p.at(this._at, this._point) && this._point.kind === 'data' ? p.strips.block[this._point.strip] : -1;
    if (this._sum > 0) {
      // out: every block again, from an opened chromosome
      if (block >= 0 && include) this._cut(section);
    } else if (block < 0) {
      // the hole: said while the ring shows the whole genome
      if (!include) v.notify(Hole);
    } else {
      const name = v.spec.model.value.systems.get(s.blocks[block].axis.system)?.name;
      const shown = s.blocks.length === 1 && include;
      if (name && !shown) this._cut(section, name);
    }
    this._done = true;
    return true;
  }

  /** `section.include` of one block (none: every block), one patch and one step of the view's history. */
  private _cut(section: string, name?: string): void {
    const v = this.viewer;
    const before = v.view.snapshot({ includes: true });
    const done = v.execute('section.include', name ? { address: section, include: [name] } : { address: section });
    if (typeof done === 'string') v.view.record(before);
    else void done.then(() => v.view.record(before));
  }
}

/** Whether the section's blocks are a chosen set (its `include`), not all of them. */
function includeOf(v: SeqViewer, section: string): boolean {
  return resolveAddress(v.spec.root.value, section.split('/').filter(Boolean))?.params?.include !== undefined;
}
