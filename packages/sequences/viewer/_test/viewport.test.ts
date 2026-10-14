/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Viewports (VC3; seqspec:interaction#viewports, seqspec:interaction#linked-sections, seqspec:interaction#focus):
// survival across document changes under both OnAxisChange values,
// a section's objects kept across a removal and a return, scroll-group joins, focus as the hull of the projection, and
// VI2's cases: one pending intent until a width, the rowZoom anchor, the animation's last frame, the `minCell` rule.
// Ranges in titles are block coordinate system indices.

import { Extents } from '@molstar/seq-core/model';
import { describe, expect, it } from 'vitest';
import {
  type FocusChange,
  manualScheduler,
  type SeqIssue,
  SeqViewerConfig,
  type SeqViewerContext,
} from '../src/index.js';
import { type Viewport, Viewports } from '../src/viewport.js';
import {
  axisOf,
  blocksDoc,
  center,
  open,
  rangedDoc,
  resolved,
  s3Doc,
  set,
  setPiece,
  vp,
  WIDTH,
} from './axis-fixtures.js';

const MODES = ['keepCenter', 'fit'] as const;
const versionOf = (v: Viewport) => v.live.version.value;
const subscribers = (s: object) => (s as { subscriberCount: number }).subscriberCount;
const view = (v: Viewport) => {
  const s = v.live.current;
  return [s.cell, s.scrollX, s.visLo[0], s.visHi[0]];
};
const issuesOf = (ctx: SeqViewerContext) => {
  const out: SeqIssue[] = [];
  ctx.events.issue.subscribe((i) => out.push(i));
  return out;
};
/** Extents of one unit's [start, stop), by its key. */
function unitExtents(ctx: SeqViewerContext, key: string, start: number, stop: number) {
  const set = ctx.graph.assembled.value?.set;
  const u = set ? set.units.key.findIndex((k) => set.strings.get(k) === key) : -1;
  const e = Extents.create();
  Extents.add(e, { kind: 1, id: u, row: -1 }, start, stop);
  return e;
}

async function opened(mode: (typeof MODES)[number]) {
  const o = await open();
  o.ctx.config.set(SeqViewerConfig.Viewport.OnAxisChange, mode);
  return { ...o, main: vp(o.ctx, 'v.main') };
}

describe('viewport survival, under both OnAxisChange values', () => {
  for (const mode of MODES) {
    describe(mode, () => {
      it('a new section opens on its focus; one without a focus fits', async () => {
        const { ctx, main } = await opened(mode);
        // Focus a's [2, 12): 10 positions in 160 px, centered.
        expect(view(main)).toEqual([16, 32, 2, 12]);
        expect(main.committed.value).toMatchObject({ origin: 'graph', blocks: [{ start: 2, stop: 12 }] });
        ctx.views.setCurrent('g');
        await ctx.graph.update();
        const one = vp(ctx, 'g.one').live.current;
        expect([one.cell, one.minCell, one.scrollX, one.visHi[0]]).toEqual([WIDTH / 20, WIDTH / 20, 0, 20]);
      });

      it("an equal axis keeps the viewport (a layer's color, a track moved)", async () => {
        const { ctx, main } = await opened(mode);
        main.zoomTo(0, 4, 8, { origin: 'api' });
        const [version, committed] = [versionOf(main), main.committed.value];
        await ctx.applyPatch(set('v.dom.l', { color: { theme: 'uniform', params: { color: '#ff0000' } } }));
        await ctx.applyPatch({ format: 'seq_view_spec_patch', ops: [{ op: 'move', ref: 'v.dom', parent: 'v.main' }] });
        expect(versionOf(main)).toBe(version);
        expect(main.committed.value).toBe(committed);
      });

      it('a changed axis keeps the center and the zoom when the center is still shown, else fits', async () => {
        const { ctx, main } = await opened(mode);
        main.zoomTo(0, 10, 14, { origin: 'api' });
        expect([main.live.current.cell, center(main)]).toEqual([28, 12]);
        // a shows 0–3 and 10–15; b, reached through the mapping, its 0–5: 16 positions and one gap.
        await ctx.applyPatch(set('v.main', { positions: 'dom' }));
        const s = main.live.current;
        const fit = (WIDTH - 28) / 16;
        expect(main.axis?.shown).toEqual(Int32Array.from([10, 6]));
        if (mode === 'keepCenter') expect([s.cell, center(main)]).toEqual([28, 12]);
        else expect([s.cell, s.minCell, s.scrollX]).toEqual([fit, fit, 0]);
        expect(main.committed.value.origin).toBe('graph');

        // The center (a's 7) is filtered out: both fit.
        await ctx.applyPatch(set('v.main', { positions: null }));
        main.zoomTo(0, 6, 8, { origin: 'api' });
        expect(center(main)).toBe(7);
        await ctx.applyPatch(set('v.main', { positions: 'dom' }));
        expect([s.cell, s.scrollX]).toEqual([fit, 0]);

        // A view of the whole filtered axis (fitted) fits the whole axis back, though its center is still shown.
        await ctx.applyPatch(set('v.main', { positions: null }));
        expect([s.cell, s.cell === s.minCell, s.scrollX]).toEqual([s.minCell, true, 0]);
        // So does one zoomed out to the whole axis by hand.
        await ctx.applyPatch(set('v.main', { positions: 'dom' }));
        main.setView(28, 0, 'api');
        main.setView(0, 0, 'api');
        expect(s.cell).toBe(s.minCell);
        await ctx.applyPatch(set('v.main', { positions: null }));
        expect([s.cell, s.scrollX]).toEqual([s.minCell, 0]);
      });

      it('a focus change is an event for the focusViewport behavior; the manager leaves the viewport', async () => {
        const { ctx, main } = await opened(mode);
        const changes: FocusChange[] = [];
        ctx.focus.changed.subscribe((c) => changes.push({ ...c }));
        const version = versionOf(main);
        await ctx.applyPatch(set('v.f', { start: 15, stop: 18 }));
        expect(changes).toMatchObject([
          { section: 'v.main', entry: { block: 0, start: 15, stop: 18 }, transitionMs: -1, origin: 'graph' },
        ]);
        expect(ctx.focus.current.value.get('v.main')).toMatchObject({ start: 15, stop: 18 });
        expect(versionOf(main)).toBe(version);
      });
    });
  }

  it('a section that leaves every view drops its states', async () => {
    const { ctx } = await open();
    const b = vp(ctx, 'g.b');
    let completed = false;
    b.live.version.subscribe({ complete: () => (completed = true) });
    await ctx.applyPatch({ format: 'seq_view_spec_patch', ops: [{ op: 'remove', ref: 'g.b' }] });
    expect(ctx.viewports.sections.value).toEqual(['v.main', 'g.one', 'g.two']);
    expect([ctx.viewports.get('g.b'), completed]).toEqual([null, true]);
  });

  it('a section that leaves and returns within one update keeps its objects; `sections` follows the sync', async () => {
    const vps = new Viewports({ scheduler: manualScheduler(), config: configOf(), active: () => true });
    const doc = s3Doc();
    vps.sync(doc);
    vps.publish();
    const main = vps.get('v.main');
    let completed = false;
    main?.live.version.subscribe({ complete: () => (completed = true) });
    // `setDocument` after `setData`: the document reset, then the new one, in one update.
    vps.sync(null);
    expect(vps.get('v.main')).toBe(null);
    vps.sync(doc);
    vps.publish();
    expect(vps.get('v.main')).toBe(main);
    await Promise.resolve();
    expect(completed).toBe(false);
    // Leaving for good: closed at the end of the update.
    vps.sync(null);
    vps.publish();
    expect(completed).toBe(false);
    await Promise.resolve();
    expect([completed, vps.sections.value]).toEqual([true, []]);
  });

  it("`sections` is emitted after every section's objects exist (through the context's fan-out)", async () => {
    const { ctx } = await open();
    const seen: boolean[] = [];
    ctx.viewports.sections.subscribe((refs) => seen.push(refs.every((r) => ctx.viewports.get(r) !== null)));
    await ctx.applyPatch({ format: 'seq_view_spec_patch', ops: [{ op: 'remove', ref: 'g.b' }] });
    expect(seen).toEqual([true, true]);
  });
});

function configOf() {
  const values = new Map<unknown, unknown>();
  return {
    get: <T>(i: { defaultValue: T }) => (values.has(i) ? (values.get(i) as T) : i.defaultValue),
    set: (i: unknown, v: unknown) => void values.set(i, v),
  };
}

describe('scroll groups', () => {
  async function group() {
    const o = await open();
    o.ctx.views.setCurrent('g');
    await o.ctx.graph.update();
    return o.ctx;
  }
  const members = (ctx: SeqViewerContext) => (ctx.scrollGroups.groups.value.get('s') ?? []).map((v) => v.ref);

  it("joins across two contexts are left on either's dispose, in both directions", async () => {
    const [a, b, c] = [await group(), await group(), await group()];
    expect(members(a)).toEqual(['g.one', 'g.two']);
    expect(a.scrollGroups.groupOf('g.two')).toBe('s');
    const joined = vp(b, 'g.one');
    const before = subscribers(joined.live.version);
    a.scrollGroups.join('s', joined);
    expect(members(a)).toEqual(['g.one', 'g.two', 'g.one']);
    b.dispose(); // the other context's dispose leaves the join here
    expect(members(a)).toEqual(['g.one', 'g.two']);

    const other = vp(c, 'g.two');
    const held = subscribers(other.live.version);
    const leave = a.scrollGroups.join('s', other);
    expect(subscribers(other.live.version)).toBe(held + 2); // the join's own and the behavior's
    a.dispose(); // this context's dispose leaves it there
    expect(subscribers(other.live.version)).toBe(held);
    expect(subscribers(joined.live.version)).toBeLessThanOrEqual(before);
    leave(); // idempotent
  });

  it('a viewport with another axis cannot join', async () => {
    const [a, b] = [await group(), await group()];
    const issues = issuesOf(a);
    a.scrollGroups.join('s', vp(b, 'g.b'));
    expect(members(a)).toEqual(['g.one', 'g.two']);
    expect(issues.map((i) => i.code)).toEqual(['svs.section.scroll_group']);
  });
});

describe('focus: the hull of the projection on the first block it lights', () => {
  it('fromExtents: every section that shows them, through the mapping too; none: false', async () => {
    const { ctx } = await open();
    ctx.views.setCurrent('g');
    await ctx.graph.update();
    const changes: FocusChange[] = [];
    ctx.focus.changed.subscribe((c) => changes.push({ ...c }));
    // a's [8, 14): on a in g.one and g.two; through the mapping, b's [0, 4) in g.b.
    expect(ctx.focus.fromExtents(unitExtents(ctx, 'a/unit', 8, 14), { origin: 'api' })).toBe(true);
    expect(changes.map((c) => [c.section, c.entry?.block, c.entry?.start, c.entry?.stop])).toEqual([
      ['g.one', 0, 8, 14],
      ['g.two', 0, 8, 14],
      ['g.b', 0, 0, 4],
    ]);
    expect(ctx.focus.fromExtents(unitExtents(ctx, 'a/unit', 0, 4), { section: 'g.b' })).toBe(false);
  });

  it("focusExtents moves a viewport to the hull on its first lit block; direct wins (b lights only b's block)", async () => {
    const { ctx } = await open();
    const m = vp(ctx, 'v.main');
    expect(m.focusExtents(unitExtents(ctx, 'a/unit', 12, 14))).toBe(true);
    expect(m.committed.value.blocks).toEqual([{ ref: 'a/unit', start: 10, stop: 16 }]);
    // b's [2, 4) is a's [12, 14) through the mapping, but the section shows b: only b's block is lit (VW5).
    expect(m.focusExtents(unitExtents(ctx, 'b/unit', 2, 4))).toBe(true);
    expect(m.committed.value.blocks).toEqual([{ ref: 'b/unit', start: 0, stop: 6 }]);
  });
});

describe('VI2: viewport deferred intents and limits', () => {
  it("a section's focus waits for a width; sized after the graph settles, the focused domain fits", async () => {
    const { ctx, size } = await open(s3Doc(), { late: true });
    const main = vp(ctx, 'v.main');
    expect(main.live.current.width).toBe(0);
    size();
    expect(view(main)).toEqual([16, 32, 2, 12]);
    expect(main.committed.value).toMatchObject({ origin: 'graph', blocks: [{ start: 2, stop: 12 }] });
  });

  it('a focus change entered on a zero-width section, moved as focusViewport does, then sized: it fits', async () => {
    const { ctx, size } = await open(s3Doc(), { late: true });
    const main = vp(ctx, 'v.main');
    ctx.focus.changed.subscribe((c) => {
      if (c.section === 'v.main' && c.entry) main.zoomTo(c.entry.block, c.entry.start, c.entry.stop, { animate: true });
    });
    await ctx.applyPatch(set('v.f', { start: 10, stop: 14 }));
    size();
    expect(view(main)).toEqual([28, 10 * 28 + 2 * 28 - WIDTH / 2, 9, 15]);
  });

  it('focusExtents before sizing wins over the pending initial focus; a pan before sizing drops the intent', async () => {
    const a = await open(s3Doc(), { late: true });
    const main = vp(a.ctx, 'v.main');
    expect(main.focusExtents(unitExtents(a.ctx, 'a/unit', 0, 20))).toBe(true);
    a.size();
    expect(view(main)).toEqual([WIDTH / 20, 0, 0, 20]);

    const b = await open(s3Doc(), { late: true });
    const other = vp(b.ctx, 'v.main');
    other.panBy(5, 0, 'pointer');
    b.size();
    expect(view(other)[0]).not.toBe(16);
    expect(other.committed.value.blocks[0]).not.toMatchObject({ start: 2, stop: 12 });
  });

  it('a hidden then shown section keeps its zoom', async () => {
    const { ctx } = await open();
    const main = vp(ctx, 'v.main');
    main.zoomTo(0, 4, 9);
    const before = view(main);
    main.setSize(0, 100, 1);
    expect(view(main)).toEqual(before);
    main.setSize(WIDTH, 100, 1);
    expect(view(main)).toEqual(before);
  });

  it("a row zoom's anchor near the bottom holds once the layout gives the rows' new height", async () => {
    const { ctx } = await open();
    const main = vp(ctx, 'v.main');
    // The viewport alone: the layout (built here) stays out of it, so its rows' heights are the test's.
    ctx.layout.get('v.main');
    if (main.axis) main.axis.ensure = null;
    main.setRows(Float64Array.of(1000, 1, 0));
    main.scrollTo(0, 900);
    main.setRowScale(2, 50);
    // Clamped against the old height until the layout sets the new one.
    expect(main.live.current.scrollY).toBe(900);
    main.setRows(Float64Array.of(2000, 1, 0));
    expect(main.live.current.scrollY).toBe((900 + 50) * 2 - 50);
    // A later move drops an anchor the layout did not take.
    main.setRowScale(1, 50);
    main.panBy(0, -10);
    main.setRows(Float64Array.of(1000, 1, 0));
    expect(main.live.current.scrollY).toBe(890);
  });

  it("an animated move (the navigator's double click) interacts until its last frame, which clears it", async () => {
    const { ctx, scheduler } = await open();
    const main = vp(ctx, 'v.main');
    main.fit();
    main.zoomCentered(28, 80, { animate: true, origin: 'behavior', durationMs: 100 });
    expect([main.live.current.animating, main.live.current.interacting]).toEqual([true, true]);
    scheduler.advance(50);
    expect(main.live.current.interacting).toBe(true);
    const commits: string[] = [];
    main.committed.subscribe((c) => commits.push(c.origin));
    scheduler.advance(60);
    expect([main.live.current.animating, main.live.current.interacting, main.live.current.cell]).toEqual([
      false,
      false,
      28,
    ]);
    expect(commits.at(-1)).toBe('behavior');
  });

  it('an animated zoom is about one view x: every position moves one way, none sweeps past and back', async () => {
    const { ctx, scheduler } = await open();
    const main = vp(ctx, 'v.main');
    main.fit();
    const s = main.live.current;
    const [c0, x0] = [s.cell, s.scrollX];
    // The full view to 28 px a position, centered at view x 150: the cell grows more than e-fold (6.4).
    main.zoomCentered(28, 150, { animate: true, origin: 'behavior', durationMs: 160 });
    const frames: [number, number][] = [[c0, x0]];
    for (let k = 0; k < 12; k++) {
      scheduler.advance(16);
      frames.push([s.cell, s.scrollX]);
    }
    const [c1, x1] = frames[frames.length - 1];
    expect(s.animating).toBe(false);
    expect(c1 / c0).toBeGreaterThan(Math.E);
    /**
     * Whether every position's view x moves one way over the frames (scroll `xs`), but for the scroll's rounding to a
     * device px (a position near the fixed x moves under 1 px a frame either way).
     */
    const oneWay = (xs: number[]) => {
      for (let p = 0; p <= 40; p++) {
        let sign = 0;
        for (let f = 1; f < frames.length; f++) {
          const d = p * frames[f][0] - xs[f] - (p * frames[f - 1][0] - xs[f - 1]);
          if (Math.abs(d) < 1) continue;
          if (sign !== 0 && Math.sign(d) !== sign) return false;
          sign = Math.sign(d);
        }
      }
      return true;
    };
    expect(oneWay(frames.map((f) => f[1]))).toBe(true);
    // Positive control: the scroll linear in the eased time, as before, sweeps positions past and back.
    const eased = frames.map(([c]) => (Math.log(c) - Math.log(c0)) / (Math.log(c1) - Math.log(c0)));
    expect(oneWay(eased.map((e) => x0 + (x1 - x0) * e))).toBe(false);
  });

  it('24 blocks in 150 px can zoom out: a quarter of the area stays for the positions, and the fit scrolls', async () => {
    const { ctx } = await open(blocksDoc(24), { late: true });
    const main = vp(ctx, ctx.viewports.sections.value[0]);
    main.setSize(150, 100, 1);
    const s = main.live.current;
    expect(main.axis?.gap).toBe(9);
    expect(s.minCell).toBeCloseTo(150 / 4 / 240, 9);
    expect(s.minCell).toBeLessThan(s.maxCell);
    expect(s.cell).toBe(s.minCell);
    expect(s.contentWidth).toBeGreaterThan(150);
    main.zoomAt(2, 75);
    expect(s.cell).toBeCloseTo((2 * 150) / 4 / 240, 9);
  });
});

describe('pieces through the viewport (reversed and cyclic pieces as AX1 and CY1 give them)', () => {
  it("a ranged unit's focus is its domain [5, 15): every position fits; committed in block indices", async () => {
    const { ctx } = await open(rangedDoc(), { late: true });
    const v = vp(ctx, ctx.viewports.sections.value[0]);
    v.setSize(WIDTH, 100, 1);
    expect(view(v)).toEqual([16, 0, 0, 10]);
    expect(v.committed.value.blocks).toEqual([{ ref: 'u', start: 5, stop: 15 }]);
    // Block [7, 9) is shown [2, 4): centered at 28 px, shown [0, 6) in view, block [5, 11).
    v.zoomTo(0, 7, 9);
    expect(view(v)).toEqual([28, 4, 0, 6]);
    expect(v.committed.value.blocks).toEqual([{ ref: 'u', start: 5, stop: 11 }]);
  });

  async function standalone(edit: (v: Parameters<typeof setPiece>[0]) => void) {
    const r = await resolved(rangedDoc(), 'w');
    edit(r.view);
    const vps = new Viewports({ scheduler: manualScheduler(), config: configOf(), active: () => true });
    const doc = rangedDoc();
    vps.sync(doc);
    vps.publish();
    const ref = r.set.strings.get(r.view.sections.ref[0]);
    vps.get(ref)?.setSize(WIDTH, 100, 1);
    vps.update(r.view, r.set, (i) => axisOf(r, i));
    return vps.get(ref) as Viewport;
  }

  it('a reversed piece: a block range zooms to its shown positions; committed as the block indices in view', async () => {
    const v = await standalone((view) => setPiece(view, 0, 0, 20, -1, 20));
    v.zoomTo(0, 15, 20);
    expect(view(v)).toEqual([28, 0, 0, 6]);
    expect(v.committed.value.blocks).toEqual([{ ref: 's/unit', start: 14, stop: 20 }]);
  });

  it('a cyclic origin: the whole piece fits; committed as its two runs', async () => {
    const v = await standalone((view) => setPiece(view, 0, 15, 5, 1, 20));
    v.zoomTo(0, 0, 20);
    expect(view(v)).toEqual([16, 0, 0, 10]);
    expect(v.committed.value.blocks).toEqual([
      { ref: 's/unit', start: 15, stop: 20 },
      { ref: 's/unit', start: 0, stop: 5 },
    ]);
  });
});
