/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The `dotPlotSelect` behavior (seqspec:interaction#selection): on a dot plot a drag selects a rectangle, its two
// ranges, the positions it spans across and down. While it lasts the pending selection holds both (its x range, and
// `plotLayer`, `yBlock`, `yStartShown`, `yEndShown`), each end following the point clamped to the plot and to the block
// the drag started on; on release the ranges, closed, replace the selection (Shift: added). Allocation-free per move.

import type { Closure } from '@molstar/seq-core/format';
import { Loci } from '@molstar/seq-core/model';
import type { SeqViewerContext } from './context.js';
import { defineBehavior } from './default-behaviors.js';
import { dotPlotAt, firstDotPlot, shownAtX } from './draw/dot-plot.js';
import { type Hit, HitRegion, HoverKind } from './hit.js';
import { Binding, DragPhase } from './input.js';
import type { SeqLociAny } from './lineage.js';
import { ofRuns } from './projection.js';

const DATA = HitRegion.data;
const START = DragPhase.start;
const END = DragPhase.end;
const CANCEL = DragPhase.cancel;
const SEQ = HoverKind.seq;

/** Whether the hit is on a `dot_plot` layer's data (edge rate: a drag's start). */
function onPlot(ctx: SeqViewerContext, hit: Hit): boolean {
  const view = ctx.graph.assembled.value?.view;
  if (!view || hit.region !== DATA || hit.layer < 0 || hit.block < 0 || hit.shown < 0) return false;
  const type = view.resolved.layers.type[hit.layer];
  return type >= 0 && view.set.strings.get(type) === 'dot_plot';
}

/** The track's first dot plot of layer k (`firstDotPlot`): the one whose overlay draws the rectangle. Edge rate. */
function plotLayerOf(ctx: SeqViewerContext, k: number): number {
  const view = ctx.graph.assembled.value?.view;
  if (!view) return k;
  const L = view.resolved.layers;
  const name = (j: number) => (L.type[j] >= 0 ? view.set.strings.get(L.type[j]) : '');
  return firstDotPlot(L.count, (j) => L.track[j], name, k);
}

export const DotPlotSelectBehavior = defineBehavior(
  'dotPlotSelect',
  'interaction',
  60,
  {
    /** The drag that selects a rectangle: a plain one. */
    dragOn: [{ button: 'primary' }] as Binding,
    /** The drag that adds its rectangle to the selection. */
    add: [{ button: 'primary', shift: true }] as Binding,
    /** What the ranges are closed by (seqspec:interaction#channels), as a click's. */
    closure: 'containing' as Closure,
  },
  'On a dot plot a drag selects a rectangle: the positions it spans across and down, as two ranges (Shift adds them).',
  ({ ctx, h, on }) => {
    const pd = ctx.selection.pending.current;
    // The pending selection is the selection manager's live record: its owner bumps it (as `select` does).
    const bump = (ctx.selection.pending as { bump?: () => void }).bump as () => void;
    /** A stacked section's lane block (the press's): the plot's both axes; -1: side by side. */
    let only = -1;
    /** The plot's entry at the press: its top and height (CSS px, the surface's), and whether the drag adds. */
    const frame = new Float64Array(2);
    let adding = false;
    const io = new Float64Array(2);
    const at = new Int32Array(2);
    let runs = new Int32Array(16);
    const axisOf = (section: number) => {
      const ref = ctx.layout.sectionRef(section);
      return ref === null ? null : ctx.layout.axis(ref);
    };
    /**
     * The ends at the point, each clamped to the plot and to its block. The point is read from the input's record: x
     * and y passed as arguments were boxed doubles, 18 bytes a move under V8's default tuning (2026-10-15).
     */
    const follow = (section: number, p: { readonly x: number; readonly y: number }) => {
      const v = axisOf(section);
      const ax = v?.axis;
      if (!v || !ax) return;
      io[0] = p.x;
      if (shownAtX(v, io, true, only, at))
        pd.endShown = at[0] === pd.block ? at[1] : at[0] > pd.block ? ax.shown[pd.block] - 1 : 0;
      io[1] = p.y;
      if (dotPlotAt(v, frame[0], frame[1], io, true, only, at))
        pd.yEndShown = at[0] === pd.yBlock ? at[1] : at[0] > pd.yBlock ? ax.shown[pd.yBlock] - 1 : 0;
    };
    const clear = () => {
      pd.active = false;
      pd.plotLayer = pd.yBlock = pd.yStartShown = pd.yEndShown = -1;
      bump();
    };
    /** Shown [min, max] of block b as seq loci; null: not a unit's. Edge rate. */
    const range = (section: number, b: number, a: number, z: number): SeqLociAny | null => {
      const ax = axisOf(section)?.axis;
      if (!ax || ax.unit[b] < 0) return null;
      const k0 = Math.min(a, z);
      const k1 = Math.max(a, z) + 1;
      let n = ax.shownToRuns(b, k0, k1, runs);
      if (2 * n > runs.length) {
        runs = new Int32Array(2 * n);
        n = ax.shownToRuns(b, k0, k1, runs);
      }
      const indices = ofRuns(runs, n);
      return { kind: 'seq-loci', set: ax.set, elements: [{ unit: ax.unit[b], indices }] } as unknown as SeqLociAny;
    };
    on(ctx.input.events.drag, (e) => {
      const p = h.forSection(e.start.section);
      if (!p) return;
      if (e.phase === START) {
        if (e.owner !== '' || !onPlot(ctx, e.start)) return;
        const add = Binding.match(p.add, e.button, e.input.modifiers);
        if (!add && !Binding.match(p.dragOn, e.button, e.input.modifiers)) return;
        const ref = ctx.layout.sectionRef(e.start.section);
        const l = ref === null ? null : ctx.layout.get(ref)?.current;
        const v = axisOf(e.start.section);
        const entry = e.start.entry;
        if (!l || !v || entry < 0) return;
        frame[0] = entry < l.pinnedCount ? l.entryTop[entry] : l.scrolledTop + l.entryTop[entry] - v.scrollY;
        frame[1] = l.entryHeight[entry];
        io[1] = e.start.y;
        only = v.axis?.stacked ? e.start.block : -1;
        if (!dotPlotAt(v, frame[0], frame[1], io, true, only, at)) return;
        e.claim('dotPlotSelect');
        adding = add;
        pd.active = true;
        pd.remove = false;
        pd.section = e.start.section;
        pd.kind = SEQ;
        pd.block = e.start.block;
        pd.startShown = pd.endShown = e.start.shown;
        pd.rowStart = pd.rowEnd = -1;
        pd.plotLayer = plotLayerOf(ctx, e.start.layer);
        pd.yBlock = at[0];
        pd.yStartShown = pd.yEndShown = at[1];
      } else if (e.owner !== 'dotPlotSelect') return;
      if (e.phase === CANCEL) return clear();
      follow(pd.section, e.input);
      if (e.phase !== END) return bump();
      const x = range(pd.section, pd.block, pd.startShown, pd.endShown);
      const y = range(pd.section, pd.yBlock, pd.yStartShown, pd.yEndShown);
      clear();
      const both = x && y ? (Loci.union(x as never, y as never) as SeqLociAny) : (x ?? y);
      if (!both) return;
      const closed = ctx.selection.close(both, p.closure);
      if (adding) ctx.selection.add(closed, 'pointer');
      else ctx.selection.set(closed, 'pointer');
    });
  },
  {
    dragOn: { action: 'dragSelect', description: 'Drag a rectangle on a dot plot with ${triggers}' },
    add: { action: 'dragAdd', description: 'Add a rectangle on a dot plot with ${triggers}' },
  },
);
