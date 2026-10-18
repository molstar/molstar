/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The Introduction's link between the Seq* viewer and Nightingale (support/seqstar-link.ts), the viewer real and
// Nightingale a fake side: the range mapping (1-based inclusive display range, 0-based half-open positions), hover and
// selection both ways, and the no-echo rule (a range or a hover that came from one side is never sent back to it).

import { createSeqViewSpec } from '@molstar/seq-core/builder';
import { createRegistries, KernelRegistry, register, resolveView } from '@molstar/seq-core/format';
import { Loci } from '@molstar/seq-core/model';
import { extentLoci, materialize, NightingaleSupported } from '@molstar/seq-nightingale';
import type { NightingaleEvent } from '@molstar/seq-nightingale/mount';
import { afterEach, describe, expect, it } from 'vitest';
import { type Opened, opened } from '../../viewer/_test/layout-fixtures.js';
import {
  fromDisplay,
  linkNightingale,
  NIGHTINGALE,
  type NightingaleSide,
  toDisplay,
} from '../src/_stories/support/seqstar-link.js';

describe('the range mapping', () => {
  it("is Nightingale's 1-based inclusive range against the viewer's 0-based half-open one, both ways", () => {
    expect(toDisplay({ start: 0, stop: 72 })).toEqual({ start: 1, end: 72 });
    expect(toDisplay({ start: 9, stop: 30 })).toEqual({ start: 10, end: 30 });
    expect(fromDisplay(1, 72)).toEqual({ start: 0, stop: 72 });
    expect(fromDisplay(10, 30)).toEqual({ start: 9, stop: 30 });
    for (const [a, b] of [
      [0, 1],
      [3, 40],
      [71, 72],
    ])
      expect(fromDisplay(toDisplay({ start: a, stop: b }).start, toDisplay({ start: a, stop: b }).end)).toEqual({
        start: a,
        stop: b,
      });
  });
  it('shows every position a brush touches, and nothing for an empty or unreadable range', () => {
    expect(fromDisplay(10.4, 20.2)).toEqual({ start: 9, stop: 21 });
    expect(fromDisplay(-3, 5)).toEqual({ start: 0, stop: 5 });
    expect(fromDisplay(7, 7)).toEqual({ start: 6, stop: 7 });
    expect(fromDisplay(8, 7)).toBeNull();
    expect(fromDisplay(Number.NaN, 7)).toBeNull();
  });
});

const reg = createRegistries();
register(reg, [KernelRegistry]);

function exampleDoc() {
  const b = createSeqViewSpec({ title: 'link' });
  const s = b.data.sequence({
    ref: 'protein',
    value: 'MSDNELKQLAEEGRKLVTALSPEDAIRWLHQNGFDVSKAIELLKGGDPSTVYYAVKIGNTEEAKQFLRDHGF',
  });
  const view = b.view({ ref: 'tracks', label: 'Tracks' }).block({ sequence: s });
  view.track({ label: '' }).ruler();
  view.track({ label: 'Residues' }).letters();
  return b.getDocument();
}

/** A Nightingale side that records what the link sends it and lets a test play what Nightingale does. */
function fake(t: Opened) {
  const assembled = t.ctx.graph.assembled.value;
  if (!assembled?.view) throw new Error('not assembled');
  const set = assembled.set;
  const { view } = resolveView(t.ctx.document.shown.value as never, set, 'tracks', reg, NightingaleSupported);
  materialize(view, set);
  const events = new Set<(e: NightingaleEvent) => void>();
  const ranges = new Set<(a: number, b: number) => void>();
  const sent = { highlights: [] as (Loci | null)[], ranges: [] as [number, number][] };
  const side: NightingaleSide = {
    view,
    set,
    highlight: (l) => sent.highlights.push(l),
    setRange: (a, b) => sent.ranges.push([a, b]),
    onRange: (f) => {
      ranges.add(f);
      return () => ranges.delete(f);
    },
    onEvent: (f) => {
      events.add(f);
      return () => events.delete(f);
    },
  };
  const event = (kind: NightingaleEvent['kind'], start = -1, stop = -1) => {
    const e: NightingaleEvent = {
      kind,
      section: 0,
      block: 0,
      index: start,
      layer: -1,
      annotation: -1,
      item: -1,
      extent: kind === 'leave' ? [] : [{ block: 0, start, stop }],
      hint: null,
    };
    for (const f of events) f(e);
  };
  return {
    side,
    sent,
    set,
    view,
    event,
    /** The user zooms or pans Nightingale to the display range. */
    zoom: (a: number, b: number) => {
      for (const f of ranges) f(a, b);
    },
  };
}

let unlink: (() => void) | null = null;
afterEach(() => {
  unlink?.();
  unlink = null;
});
async function linked() {
  const t = await opened(exampleDoc(), { presets: false });
  await t.ctx.graph.update();
  const f = fake(t);
  // The timer of the deferred put-back of a selection under Nightingale's own hover: run now.
  unlink = linkNightingale(t.ctx, f.side, (run) => run());
  const section = t.ctx.viewports.sections.value[0];
  const vp = t.ctx.viewports.get(section);
  if (!vp) throw new Error('no viewport');
  return { t, f, vp, ctx: t.ctx };
}
const holds = (a: Loci | null | undefined, b: Loci) => !!a && Loci.areEqual(a, b);

describe('the hover', () => {
  it("the viewer's lights Nightingale's, and clears it", async () => {
    const { f, ctx } = await linked();
    const lit = extentLoci(f.view, f.set, [{ block: 0, start: 4, stop: 6 }]) as Loci;
    f.sent.highlights.length = 0;
    ctx.hover.setLoci(lit, 'api');
    expect(f.sent.highlights).toHaveLength(1);
    expect(holds(f.sent.highlights[0], lit)).toBe(true);
    ctx.hover.setLoci(null, 'api');
    expect(f.sent.highlights.at(-1)).toBeNull();
  });

  it("Nightingale's lights the viewer's, as its own move, and is not sent back", async () => {
    const { f, ctx } = await linked();
    f.sent.highlights.length = 0;
    f.event('hover', 10, 11);
    expect(ctx.hover.live.current.owner).toBe(NIGHTINGALE);
    expect(holds(ctx.hover.toLoci(), extentLoci(f.view, f.set, [{ block: 0, start: 10, stop: 11 }]) as Loci)).toBe(
      true,
    );
    expect(f.sent.highlights).toEqual([]);
    f.event('leave');
    expect(ctx.hover.live.current.owner).toBe('');
    expect(f.sent.highlights).toEqual([]);
  });
});

describe('the selection', () => {
  it("Nightingale's click selects in the viewer, the same click again clears it, and the viewer's selection lights Nightingale", async () => {
    const { f, ctx } = await linked();
    const first = extentLoci(f.view, f.set, [{ block: 0, start: 4, stop: 12 }]) as Loci;
    f.sent.highlights.length = 0;
    f.event('click', 4, 12);
    expect(ctx.selection.current.value.origin).toBe(NIGHTINGALE);
    expect(holds(ctx.selection.current.value.seq, first)).toBe(true);
    // Its click is shown in Nightingale as the viewer's selection is: its own highlight does not hold a click.
    expect(holds(f.sent.highlights.at(-1), first)).toBe(true);
    f.event('click', 20, 22);
    expect(
      holds(ctx.selection.current.value.seq, extentLoci(f.view, f.set, [{ block: 0, start: 20, stop: 22 }]) as Loci),
    ).toBe(true);
    f.event('click', 20, 22);
    expect(ctx.selection.current.value.seq).toBeNull();
    expect(f.sent.highlights.at(-1)).toBeNull();
  });

  it("a hover over a selection lights both in Nightingale, a hover of Nightingale's own puts the selection back", async () => {
    const { f, ctx } = await linked();
    const picked = extentLoci(f.view, f.set, [{ block: 0, start: 4, stop: 6 }]) as Loci;
    ctx.selection.set(picked, 'api');
    expect(holds(f.sent.highlights.at(-1), picked)).toBe(true);
    ctx.hover.setLoci(extentLoci(f.view, f.set, [{ block: 0, start: 30, stop: 31 }]) as Loci, 'api');
    const both = extentLoci(f.view, f.set, [
      { block: 0, start: 4, stop: 6 },
      { block: 0, start: 30, stop: 31 },
    ]) as Loci;
    expect(holds(f.sent.highlights.at(-1), both)).toBe(true);
    ctx.hover.setLoci(null, 'api');
    f.sent.highlights.length = 0;
    f.event('hover', 40, 41);
    expect(
      holds(
        f.sent.highlights.at(-1),
        extentLoci(f.view, f.set, [
          { block: 0, start: 40, stop: 41 },
          { block: 0, start: 4, stop: 6 },
        ]) as Loci,
      ),
    ).toBe(true);
  });
});

describe('the current view', () => {
  it("the viewer's zoom moves Nightingale to the same positions, once, in its 1-based range", async () => {
    const { f, vp } = await linked();
    f.sent.ranges.length = 0;
    vp.zoomTo(0, 9, 30, { animate: false });
    const [last] = f.sent.ranges.slice(-1);
    // The viewer shows whole cells: the range it shows covers the one asked for.
    expect(last[0]).toBeLessThanOrEqual(10);
    expect(last[1]).toBeGreaterThanOrEqual(30);
    const shown = vp.live.current;
    expect([last[0] - 1, last[1]]).toEqual([shown.visLo[0], shown.visHi[0]]);
    const n = f.sent.ranges.length;
    vp.zoomTo(0, 9, 30, { animate: false });
    expect(f.sent.ranges).toHaveLength(n);
  });

  it("Nightingale's zoom moves the viewer, and the viewer's settling there is not sent back", async () => {
    const { f, vp } = await linked();
    f.sent.ranges.length = 0;
    f.zoom(41, 60);
    const shown = vp.live.current;
    expect(shown.visLo[0]).toBeLessThanOrEqual(40);
    expect(shown.visHi[0]).toBeGreaterThanOrEqual(60);
    expect(vp.live.current.origin).toBe(NIGHTINGALE);
    expect(f.sent.ranges).toEqual([]);
    // Its echo, the range it was set to, is not a move either.
    f.zoom(41, 60);
    expect(f.sent.ranges).toEqual([]);
  });

  it("each side's next move after the other's is sent, so the two never stick", async () => {
    const { f, vp } = await linked();
    f.zoom(41, 60);
    f.sent.ranges.length = 0;
    vp.zoomTo(0, 2, 20, { animate: false });
    expect(f.sent.ranges).toHaveLength(1);
    vp.panBy(30, 0, 'pointer');
    expect(f.sent.ranges).toHaveLength(2);
    const before = vp.live.current.visLo[0];
    f.zoom(5, 25);
    expect(vp.live.current.visLo[0]).not.toBe(before);
    expect(f.sent.ranges).toHaveLength(2);
  });
});
