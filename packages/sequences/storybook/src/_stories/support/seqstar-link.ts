/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Two renderings of one document linked both ways, the Seq* viewer's and Nightingale's (the Introduction page's pair): the
// hover and the selection of one light the other, and a zoom or a pan of one moves the other to the same range. The
// shape of the Mol* adapter's `linkViews` (Prism's) (seqspec:libraries#the-mol-adapter): Nightingale is a host of the viewer, its
// events are the viewer's `host:nightingale` moves, and nothing a side hears from the other is sent back to it.
//
// The viewer's replays and carries are no one's move (seqspec:interaction#replays-and-carries): the viewport's and the
// selection's replay on subscribing is skipped, and a range is sent only where it differs from the last one sent either
// way (so Nightingale's echo of a range, and the viewer's settling on the range it was given, send nothing). One
// section, its first block: a range is that block's positions.

import type { ResolvedView } from '@molstar/seq-core/format';
import { Loci, type SeqSet } from '@molstar/seq-core/model';
import { extentLoci } from '@molstar/seq-nightingale';
import type { NightingaleEvent } from '@molstar/seq-nightingale/mount';
import type { SeqViewerContext } from '@molstar/seq-viewer';

/** The origin Nightingale's moves carry in the viewer. */
export const NIGHTINGALE = 'host:nightingale';

/** What the link needs of Nightingale's side: `mount`'s, over the first section's manager (a fake in tests). */
export interface NightingaleSide {
  /** The view and set it draws: its events' extents are loci of this set, which is the viewer's. */
  readonly view: ResolvedView;
  readonly set: SeqSet;
  /** `Mounted.highlight`: lights the loci (null: clears) in every manager. */
  highlight(loci: Loci | null): void;
  /** Sets the manager's display range, 1-based and inclusive, as a section's zoom does (`display-start`, `display-end`). */
  setRange(start: number, end: number): void;
  /** Nightingale's own zoom or pan, a user's: its `display-start` and `display-end`; never the one `setRange` sets. */
  onRange(listener: (start: number, end: number) => void): () => void;
  /** `mount`'s events: hover, click and leave, in the view's terms. */
  onEvent(listener: (e: NightingaleEvent) => void): () => void;
}

/** Positions [start, stop) of the block, 0-based as the viewer's. */
export interface Range {
  readonly start: number;
  readonly stop: number;
}
const sameRange = (a: Range | null, b: Range | null) =>
  a === b || (!!a && !!b && a.start === b.start && a.stop === b.stop);

/** Nightingale's 1-based inclusive display range as the block's positions; null: nothing to show (empty or not numbers). */
export function fromDisplay(start: number, end: number): Range | null {
  if (!Number.isFinite(start) || !Number.isFinite(end)) return null;
  // A brush ends inside a position: every position it touches is shown.
  const range = { start: Math.max(0, Math.floor(start) - 1), stop: Math.ceil(end) };
  return range.stop > range.start ? range : null;
}
/** The block's positions as Nightingale's display range. */
export const toDisplay = (r: Range) => ({ start: r.start + 1, end: r.stop });

/** The part of a selection Nightingale shows: its highlight is one loci, the first part the selection holds. */
function selected(ctx: SeqViewerContext): Loci | null {
  const s = ctx.selection.current.value;
  return (s.seq ?? s.items ?? s.aligned ?? s.columns ?? s.rows ?? s.boundaries) as Loci | null;
}
/** Both, when they are of one kind (a union); else the first (a hover over a selection of another kind). */
function together(a: Loci | null, b: Loci | null): Loci | null {
  if (!a || !b) return a ?? b;
  return a.kind === b.kind && a.set === b.set && !('alignment' in a && 'alignment' in b && a.alignment !== b.alignment)
    ? Loci.union(a, b)
    : a;
}

/**
 * Links `ctx` and `side` (both ways), until the returned unlink. The viewer's hover and selection are shown in
 * Nightingale's highlight; Nightingale's hover and click are the viewer's hover and selection (a click on what its own
 * click selected clears it); the two zooms follow each other.
 */
export function linkNightingale(
  ctx: SeqViewerContext,
  side: NightingaleSide,
  defer: (run: () => void) => void = (f) => setTimeout(f),
) {
  const stops: (() => void)[] = [];
  const section = () => ctx.viewports.sections.value[0];

  // The viewer's hover and selection in Nightingale. Nightingale shows its own hover on its own, and a selection the
  // viewer holds is put back after it: its hover is not sent again, but the union with the selection is.
  let told = false;
  const show = () => {
    const own = told || ctx.hover.live.current.owner === NIGHTINGALE;
    const sel = selected(ctx);
    if (own && !sel) return;
    side.highlight(together(ctx.hover.toLoci() as Loci | null, sel));
  };
  const hover = ctx.hover.live.version.subscribe(show);
  const selection = ctx.selection.current.subscribe(show);
  stops.push(
    () => hover.unsubscribe(),
    () => selection.unsubscribe(),
  );

  // Nightingale's hover and click in the viewer.
  stops.push(
    side.onEvent((e) => {
      const loci = e.kind === 'leave' ? null : extentLoci(side.view, side.set, e.extent);
      if (e.kind === 'click') {
        const held = ctx.selection.current.value;
        if (!loci) return;
        if (held.origin === NIGHTINGALE && held.seq && Loci.areEqual(held.seq, loci)) ctx.selection.clear(NIGHTINGALE);
        else ctx.selection.set(loci, NIGHTINGALE);
        return;
      }
      told = true;
      ctx.hover.setLoci(loci, NIGHTINGALE);
      told = false;
      // Nightingale's own highlight follows its events and takes the selection with it: put it back after.
      if (selected(ctx)) defer(show);
    }),
  );

  // The two zooms. `last`: the range last sent either way.
  let last: Range | null = null;
  let replay = true;
  const move = ctx.viewports.get(section())?.live.version.subscribe(() => {
    const first = replay;
    replay = false;
    const live = ctx.viewports.get(section())?.live.current;
    if (first || !live || live.origin === NIGHTINGALE || live.visHi.length === 0) return;
    const range = { start: live.visLo[0], stop: live.visHi[0] };
    if (range.stop <= range.start || sameRange(range, last)) return;
    last = range;
    const d = toDisplay(range);
    side.setRange(d.start, d.end);
  });
  if (move) stops.push(() => move.unsubscribe());
  stops.push(
    side.onRange((start, end) => {
      const range = fromDisplay(start, end);
      if (!range || sameRange(range, last)) return;
      last = range;
      ctx.viewports.get(section())?.zoomTo(0, range.start, range.stop, { origin: NIGHTINGALE, animate: false });
    }),
  );
  return () => {
    for (const stop of stops) stop();
    if (ctx.disposed) return;
    // Nothing of ours stays lit in either.
    ctx.hover.clear(NIGHTINGALE);
    side.highlight(null);
  };
}
