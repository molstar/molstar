/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { resolveAddress } from '@molstar/core/repr/tree';
import type { ModifierState } from '@molstar/seqx-viewer/bindings';
import type { HitTarget } from '@molstar/seqx-viewer/models';
import type { SeqViewer } from './viewer.js';

// The zoom of a ring nobody follows is a cut (_spec/README.md "Input", The ring's zoom is a cut): a ring shows its whole axis, so its window is idle; the zoom gesture opens the chromosome under the
// pointer (`section.include` of its block: in) or shows every block again (out), one cut a wheel gesture (the deltas
// past `navigate.cutNotch` px, then nothing until the wheel rests `navigate.cutRestMs`). A zoom that changes nothing says
// once a gesture what the ring shows; outside the ring's box it is the bindings'. A ring a section follows keeps the zoom of its window. A cut is the command the stem
// menu's "Show only" runs, one patch and one step of the view's history: nothing is drawn or kept per frame.

const Whole = 'The ring shows the whole genome: open a chromosome from its label’s menu or the square';
const Filtered = 'The ring shows the chromosomes its view selects: open one from its label’s menu or the square';
/** What an opened ring shows, by its blocks' count and its first's name. */
const shows = (n: number, first: string | undefined) =>
  n === 1
    ? `The ring shows ${first}: wheel out over it to show every chromosome`
    : `The ring shows ${n} chromosomes: wheel out over one to show every chromosome`;

export class RingCut {
  /** When the last wheel of the gesture came (ms); the gesture is over after a rest. */
  private _last = Number.NEGATIVE_INFINITY;
  /** The gesture's deltas summed (negative: in). */
  private _sum = 0;
  /** The gesture has cut or said: nothing more until it rests. */
  private _done = false;

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
    const sections = v.view.resolved.value?.sections;
    const s = sections?.find((x) => x.address === section);
    const p = v.placementOf(section);
    // a ring idle as the view model has it: following none and followed by none (a detached follower's window moves)
    if (!s?.blocks || !p?.warpOf || !p.blockAt || s.follow || sections!.some((x) => x.follow?.leader === section))
      return false;
    // outside the ring's box (its margin) the zoom is the bindings' own
    const block = p.blockAt(x, y);
    if (block === undefined) return false;
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
    this._done = true;
    const name = (k: number) => v.spec.model.value.systems.get(s.blocks![k].axis.system)?.name;
    const include = includeOf(v, section);
    // an author's filter is no opened chromosome: the ring says so and keeps it
    if (include === 'filter') return this._say(Filtered);
    if (this._sum > 0) {
      // out: every block again, from an opened chromosome
      if (include === 'none') return this._say(Whole);
      if (block >= 0) this._cut(section);
      else this._say(shows(s.blocks.length, name(0)));
    } else if (block < 0) {
      this._say(include === 'none' ? Whole : shows(s.blocks.length, name(0)));
    } else {
      const to = name(block);
      if (!to) return true;
      if (include === 'list' && s.blocks.length === 1) this._say(shows(1, to));
      else this._cut(section, to);
    }
    return true;
  }

  private _say(text: string): true {
    this.viewer.notify(text);
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

/**
 * What the section's `include` is: `list` (chromosomes chosen by name: an opened ring), `filter` (a pattern or a
 * `{ min_length }`, the author's), `none` (every block).
 */
function includeOf(v: SeqViewer, section: string): 'none' | 'list' | 'filter' {
  const inc = resolveAddress(v.spec.root.value, section.split('/').filter(Boolean))?.params?.include;
  if (inc === undefined) return 'none';
  return Array.isArray(inc) && inc.every((n) => typeof n === 'string' && !/[*?[\]{}]/.test(n)) ? 'list' : 'filter';
}
