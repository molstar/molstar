/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// CG3 (seqspec:interaction#linked-sections): the ring's window on the linear section linked by `follow`: a sector of
// the linear view, following it (a zoom, a pan, across a rotated cyclic origin); dragging it pans that section, a grip
// zooms it with the other end held; the cursor's part; drawn as sectors with grips on the rim.

import { createSeqViewSpec, translate } from '@molstar/seq-core/builder';
import { describe, expect, it } from 'vitest';
import { Dirty, Pass, SectionComposer, SkinToken } from '../src/draw/index.js';
import { placedRingLabels } from '../src/draw/ring-labels.js';
import { Hit, HitRegion, HoverKind } from '../src/hit.js';
import { KernelRenderers, type LayerRenderer } from '../src/index.js';
import { DefaultPresets } from '../src/presets.js';
import { Ring } from '../src/ring.js';
import { RING_WINDOW, ringBounds, ringWindowAt, ringWindowFrame } from '../src/ring-window.js';
import { Op, type RecordingBackend } from '../src/testing/index.js';
import type { Viewport } from '../src/viewport.js';
import { WindowKind, WindowPart } from '../src/windows.js';
import { recording } from './composer-fixtures.js';
import { driven } from './input-fixtures.js';
import { N, type RingWindowOptions, ringWindowDoc } from './ring-window-fixtures.js';

const registry = [{ presets: DefaultPresets.filter((p) => p.name === 'circular') }];

async function open(origin = 0, o: RingWindowOptions = {}) {
  const t = await driven(ringWindowDoc(origin, N, o), {
    registry,
    renderers: KernelRenderers as LayerRenderer[],
    size: { ring: [700, 500], lin: [1000, 200], lin2: [800, 200] },
  });
  const lin = t.ctx.viewports.get('lin') as Viewport;
  const ring = t.ctx.layout.sectionIndex('ring');
  const g = new Float64Array(RING_WINDOW);
  /** The window's ends as ring positions (the ring's shown positions, unwrapped past the block's end). */
  const ends = () => {
    expect(ringWindowFrame(t.ctx, ring, g)).toBe(true);
    const px = g[6];
    return [(g[0] - g[2]) / px, (g[1] - g[2]) / px];
  };
  /**
   * The surface point at ring position `pos` (shown, of `block`: default the last, the linked one's) and depth `y`
   * (unrolled, CSS px below the rings' top).
   */
  const at = (pos: number, depth: number, block = -1) => {
    const l = t.ctx.layout.get('ring')?.current;
    const av = t.ctx.layout.axis('ring');
    if (!l || !av?.axis) throw new Error('no ring');
    const b = block < 0 ? av.axis.count - 1 : block;
    const r = new Ring().from(l, av, 1);
    const out = new Float64Array(2);
    r.toSurface(av.blockX[b] + pos * av.span[b] * av.cell, depth, out);
    const rb = new Float64Array(2);
    ringBounds(l, rb);
    return { x: out[0], y: out[1], rim: rb[0] };
  };
  /** The linear view's first and last positions (fractional). */
  const shown = () => {
    const s = lin.live.current;
    const area = s.width - (lin.axis?.insetX ?? 0);
    return [s.scrollX / s.cell, (s.scrollX + area) / s.cell];
  };
  return { ...t, lin, ring, g, ends, at, shown };
}

describe("the ring's window", () => {
  it("is the linear section's view on the ring: the whole ring at its fit, then its zoom and pan", async () => {
    const t = await open();
    const [a0, a1] = t.ends();
    expect(a0).toBeCloseTo(0, 3);
    expect(a1).toBeCloseTo(N, 3);
    t.lin.setView(10, 1000, 'api');
    const [lo, hi] = t.shown();
    const [b0, b1] = t.ends();
    expect(b0).toBeCloseTo(lo, 3);
    expect(b1).toBeCloseTo(hi, 3);
    expect(b1 - b0).toBeLessThan(N / 2);
  });

  it('wraps across a rotated cyclic origin', async () => {
    // Shown 0 is index 550: the linear view of [530, 570) is the ring's [580, 620): past the block's end.
    const t = await open(550);
    const s = t.lin.live.current;
    const area = s.width - (t.lin.axis?.insetX ?? 0);
    t.lin.setView(area / 40, 530 * (area / 40), 'api');
    const [a, b] = t.ends();
    expect(a).toBeCloseTo(580, 1);
    expect(b).toBeCloseTo(620, 1);
    expect(t.g[4]).toBe(1);
  });

  it('pans the linear section when dragged, and a grip zooms it with the other end held', async () => {
    const t = await open();
    t.lin.setView(10, 1000, 'api');
    const [lo, hi] = t.shown();
    const cell = t.lin.live.current.cell;
    const mid = (lo + hi) / 2;
    const p = t.at(mid, 30);
    expect(ringWindowAt(t.ctx, t.ring, p, t.g)).toBe(WindowPart.pan);
    const q = t.at(mid + 40, 30);
    t.drag('ring', p.x, p.y, q.x, q.y);
    const [lo2, hi2] = t.shown();
    expect(t.lin.live.current.cell).toBeCloseTo(cell, 6);
    expect(lo2 - lo).toBeCloseTo(40, 0);
    expect(hi2 - hi).toBeCloseTo(40, 0);
    // The end grip on the rim, dragged out by 60 positions: the start held, the view wider.
    const e = t.at(hi2, t.at(0, 0).rim);
    expect(ringWindowAt(t.ctx, t.ring, e, t.g)).toBe(WindowPart.right);
    const f = t.at(hi2 + 60, t.at(0, 0).rim);
    t.drag('ring', e.x, e.y, f.x, f.y);
    const [lo3, hi3] = t.shown();
    expect(lo3).toBeCloseTo(lo2, 0);
    expect(hi3).toBeCloseTo(hi2 + 60, 0);
  });

  it('lights its part under the pointer for the cursor, and leaves a press outside it to the rings', async () => {
    const t = await open();
    t.lin.setView(10, 1000, 'api');
    const [lo, hi] = t.shown();
    const p = t.at((lo + hi) / 2, 30);
    t.pointer('move', 'ring', p.x, p.y);
    expect(t.ctx.ui.window.current.kind).toBe(WindowKind.ring);
    expect(t.ctx.ui.window.current.part).toBe(WindowPart.pan);
    expect(t.ctx.ui.cursor('ring').value).toBe('grab');
    const out = t.at((hi + 50) % N, 30);
    expect(ringWindowAt(t.ctx, t.ring, out, t.g)).toBe(WindowPart.none);
    const before = t.lin.live.current.scrollX;
    t.drag('ring', out.x, out.y, t.at((hi + 90) % N, 30).x, t.at((hi + 90) % N, 30).y);
    expect(t.lin.live.current.scrollX).toBe(before);
  });
});

/** A pass's sectors: [r0, r1, a0, a1, color, alpha]. */
function overlaySectors(r: RecordingBackend): number[][] {
  const b = r.sectorItems;
  const out: number[][] = [];
  let at = -1;
  for (let k = 0; k < r.opCount; k++) {
    const code = r.ops[4 * k];
    if (code === Op.beginPass) at = r.args[r.ops[4 * k + 3]];
    if (code !== Op.sectors || at !== Pass.overlay) continue;
    for (let i = r.ops[4 * k + 1]; i < r.ops[4 * k + 1] + r.ops[4 * k + 2]; i++)
      out.push([b.r0[i], b.r1[i], b.a0[i], b.a1[i], b.color[i], b.alpha[i]]);
  }
  return out;
}

describe("the ring's window, drawn", () => {
  it('is a `--seq-window` sector over the rings with its sides, its arc on the rim and two grips; at the fit, the grips only', async () => {
    const t = await open(550);
    const composer = new SectionComposer(t.ctx, 'ring');
    const rec = recording();
    const frame = () => {
      rec.reset();
      composer.frame(rec, Dirty.overlay | Dirty.tracks);
      return overlaySectors(rec);
    };
    const by = (all: number[][], c: number) => all.filter((x) => x[4] === c);
    let all = frame();
    expect(by(all, SkinToken.window)).toHaveLength(0);
    expect(by(all, SkinToken.grip)).toHaveLength(2);
    // [530, 570) of the linear view: the ring's [580, 620), across the rotated origin in two parts.
    const s = t.lin.live.current;
    const area = s.width - (t.lin.axis?.insetX ?? 0);
    t.lin.setView(area / 40, 530 * (area / 40), 'api');
    all = frame();
    const fill = by(all, SkinToken.window);
    expect(fill).toHaveLength(2);
    const TAU = 2 * Math.PI;
    const sweep = fill.reduce((n, x) => n + (x[3] - x[2]), 0);
    expect(sweep).toBeCloseTo((40 / N) * TAU, 1);
    // Its two sides and two grips in the grips' ink, nothing in a selection's.
    expect(by(all, SkinToken.grip)).toHaveLength(4);
    expect(by(all, SkinToken.selection).length + by(all, SkinToken.selectionBorder).length).toBe(0);
    // A lit grip: the cursor's part over the start's.
    t.ctx.ui.setWindow(t.ring, WindowKind.ring, WindowPart.left);
    all = frame();
    expect([by(all, SkinToken.grip).length, by(all, SkinToken.gripHover).length]).toEqual([3, 1]);
    composer.dispose();
  });
});

describe('labels on rings', () => {
  it('hover and select their item, as the item does', async () => {
    const t = await open();
    const composer = new SectionComposer(t.ctx, 'ring');
    composer.frame(recording(), Dirty.all);
    const layout = t.ctx.layout.get('ring')?.current as object;
    const labels = placedRingLabels.get(layout);
    expect(labels?.n).toBe(3);
    const L = labels as NonNullable<typeof labels>;
    const items = new Set<number>();
    for (let k = 0; k < L.n; k++) {
      // The text's middle: beside its anchor, on its side.
      const x = L.pts[10 * k + 8] + L.side[k] * 8;
      const y = L.pts[10 * k + 9];
      const hit = Hit.create();
      expect(t.ctx.layout.hitTest('ring', x, y, hit)).toBe(true);
      expect([hit.region, hit.kind, hit.annotation, hit.layer]).toEqual([
        HitRegion.data,
        HoverKind.item,
        0,
        L.item[3 * k],
      ]);
      items.add(hit.item);
      t.pointer('move', 'ring', x, y);
      expect(t.ctx.hover.live.current.kind).toBe(HoverKind.item);
      expect(t.ctx.labels.current.value).not.toBeNull();
    }
    expect(items.size).toBe(3);
    // A click on the last selects its item; off the labels, outside the rings, nothing.
    const x = L.pts[10 * (L.n - 1) + 8] + L.side[L.n - 1] * 8;
    t.click('ring', x, L.pts[10 * (L.n - 1) + 9]);
    expect(t.ctx.selection.current.value.items).not.toBeNull();
    const hit = Hit.create();
    t.ctx.layout.hitTest('ring', 2, 2, hit);
    expect(hit.kind).toBe(HoverKind.none);
    composer.dispose();
  });
});

/** A linear section's center element (index with its fraction) and px a shown position. */
function centerOf(vp: Viewport) {
  const out = new Float64Array(2);
  vp.centerElement(out);
  return { u: out[0], px: out[1] };
}

describe('a ring is a pass-through for its followers, and leads only through its window', () => {
  it("the linear detail opens at its fit, not at the ring's cell", async () => {
    const t = await open();
    const s = t.lin.live.current;
    expect(s.scrollX).toBe(0);
    expect(s.cell).toBeCloseTo(s.minCell, 6);
    expect(t.shown()[1]).toBeCloseTo(N, 1);
  });

  it("two followers of a ring follow each other through it; neither is snapped back to the ring's view", async () => {
    const t = await open(0, { followers: 2 });
    const l2 = t.ctx.viewports.get('lin2') as Viewport;
    t.lin.setView(10, 3000, 'api');
    const c1 = centerOf(t.lin);
    expect(centerOf(l2).u).toBeCloseTo(c1.u, 1);
    expect(centerOf(l2).px).toBeCloseTo(c1.px, 6);
    l2.setView(5, 1000, 'api');
    const c2 = centerOf(l2);
    expect(centerOf(t.lin).u).toBeCloseTo(c2.u, 1);
    expect(t.lin.live.current.cell).toBeCloseTo(5, 6);
  });

  it("the ring's + and - zoom the detail about the window's center; the ring's own viewport leads nothing", async () => {
    const t = await open();
    t.lin.setView(10, 3000, 'api');
    const [lo, hi] = t.shown();
    t.focus('ring');
    const origins: string[] = [];
    t.lin.committed.subscribe((c) => origins.push(c.origin));
    origins.length = 0;
    t.key('ring', '+');
    t.advance(2000);
    // A key's move, as the linear section's own `+` (GG-04).
    expect([t.lin.live.current.origin, ...new Set(origins)]).toEqual(['keyboard', 'keyboard']);
    origins.length = 0;
    const [lo2, hi2] = t.shown();
    expect((lo2 + hi2) / 2).toBeCloseTo((lo + hi) / 2, 0);
    expect(hi2 - lo2).toBeLessThan(hi - lo);
    t.key('ring', '-');
    t.key('ring', '-');
    const [lo3, hi3] = t.shown();
    expect(hi3 - lo3).toBeGreaterThan(hi - lo);
    expect((lo3 + hi3) / 2).toBeCloseTo((lo + hi) / 2, 0);
    // The ring's hidden viewport moved and resized: the detail stays.
    const before = [t.lin.live.current.cell, t.lin.live.current.scrollX];
    const rv = t.ctx.viewports.get('ring') as Viewport;
    rv.setView(3, 200, 'api');
    rv.setSize(640, 480, 1);
    expect([t.lin.live.current.cell, t.lin.live.current.scrollX]).toEqual(before);
  });
});

describe("the ring's window, its edges", () => {
  it('on a ring of two blocks, catches presses on its own block only', async () => {
    const t = await open(0, { before: true });
    t.lin.setView(10, 3000, 'api');
    const [lo, hi] = t.shown();
    const mid = (lo + hi) / 2;
    expect(ringWindowAt(t.ctx, t.ring, t.at(mid, 30), t.g)).toBe(WindowPart.pan);
    // The same angle along the other block: never the window.
    for (const pos of [mid, 10, 300, 590])
      expect([pos, ringWindowAt(t.ctx, t.ring, t.at(pos, 30, 0), t.g)]).toEqual([pos, WindowPart.none]);
    // A drag on the window pans by its travel from the first move (no jump of a span).
    const p = t.at(mid, 30);
    const q = t.at(mid + 20, 30);
    t.drag('ring', p.x, p.y, q.x, q.y);
    expect(t.shown()[0] - lo).toBeCloseTo(20, 0);
  });

  it("a grip drag stops at the detail's content and at its least zoom, the other end held", async () => {
    // The detail is [100, 400] of the ring's molecule, at [200, 300).
    const t = await open(0, { range: [100, 400] });
    const area = t.lin.live.current.width - (t.lin.axis?.insetX ?? 0);
    t.lin.setView(area / 100, 100 * (area / 100), 'api');
    const [a, b] = t.ends();
    expect([a, b].map(Math.round)).toEqual([200, 300]);
    const rim = t.at(0, 0).rim;
    const l = t.at(a, rim);
    const to = t.at(50, rim);
    expect(ringWindowAt(t.ctx, t.ring, l, t.g)).toBe(WindowPart.left);
    t.drag('ring', l.x, l.y, to.x, to.y);
    const [a2, b2] = t.ends();
    expect(Math.round(b2)).toBe(300);
    expect(a2).toBeGreaterThanOrEqual(99.5);
    // The end grip out past the content's end (400): it stops there, the start held.
    const r = t.at(b2, rim);
    const far = t.at(520, rim);
    t.drag('ring', r.x, r.y, far.x, far.y);
    const [a3, b3] = t.ends();
    expect(Math.round(a3)).toBe(Math.round(a2));
    expect(b3).toBeLessThanOrEqual(401.5);
  });

  it('a pan across a rotated origin moves the detail by its travel', async () => {
    const t = await open(550);
    const area = t.lin.live.current.width - (t.lin.axis?.insetX ?? 0);
    t.lin.setView(area / 40, 520 * (area / 40), 'api');
    const [lo] = t.shown();
    // The window is the ring's [570, 610): dragged by 25 it crosses the ring's origin (shown 600).
    const p = t.at(590, 30);
    const q = t.at(615, 30);
    t.drag('ring', p.x, p.y, q.x, q.y);
    expect(t.shown()[0] - lo).toBeCloseTo(25, 0);
  });
});

describe('labels on rings after a change', () => {
  it('hit nothing between a relayout and the next frame, then their items at the new places (a rotated origin)', async () => {
    const t = await open();
    const composer = new SectionComposer(t.ctx, 'ring');
    composer.frame(recording(), Dirty.all);
    const layout = t.ctx.layout.get('ring')?.current as object;
    const L = placedRingLabels.get(layout) as NonNullable<ReturnType<typeof placedRingLabels.get>>;
    const where = (k: number) => [L.pts[10 * k + 8] + L.side[k] * 8, L.pts[10 * k + 9]] as const;
    const before = Array.from({ length: L.n }, (_, k) => where(k));
    const r = await t.ctx.applyPatch({
      format: 'seq_view_spec_patch',
      ops: [{ op: 'set', ref: 'ring.b', params: { origin: 300 } }],
    } as never);
    expect(r.issues).toEqual([]);
    const hit = Hit.create();
    // Laid out again, not yet drawn: the kept labels are a frame behind and hit nothing.
    t.ctx.layout.get('ring');
    for (const [x, y] of before) {
      t.ctx.layout.hitTest('ring', x, y, hit);
      expect(hit.kind).toBe(HoverKind.none);
    }
    composer.frame(recording(), Dirty.all);
    const after = placedRingLabels.get(t.ctx.layout.get('ring')?.current as object) as typeof L;
    expect(after.n).toBe(3);
    for (let k = 0; k < after.n; k++) {
      const [x, y] = [after.pts[10 * k + 8] + after.side[k] * 8, after.pts[10 * k + 9]];
      expect(t.ctx.layout.hitTest('ring', x, y, hit)).toBe(true);
      expect([hit.kind, hit.layer]).toEqual([HoverKind.item, after.item[3 * k]]);
    }
    composer.dispose();
  });
});

describe("the ring's window through a codon step", () => {
  /** Case 14's gene on a ring and its translation as the linear detail (a route of span 3). */
  async function gene() {
    const b = createSeqViewSpec();
    const GENE =
      'GAATTCGCCACC' +
      'ATGGCTTCCAAGGAACTGACCCCACGTGACGGCC' +
      'GTAAGTCTTTCCTAACTTCTGTTTCCACAG' +
      'ATTGGTACTTCAACAGCCAGGTGGAGTAA' +
      'GCGGCCGCATAAAGT';
    const g = b.data.sequence({ ref: 'gene', alphabet: 'dna', value: GENE });
    const protein = translate(b.data, {
      ref: 'protein',
      from: g,
      cds: [
        [12, 46],
        [76, 105],
      ],
    });
    const v = b.view({ ref: 'v', label: 'V' });
    const ring = v.section({ ref: 'ring', style: 'circular' }).block({ sequence: g });
    ring.track({ label: '' }).ruler();
    v.section({ ref: 'lin', follow: ring.handle }).block({ sequence: protein }).track().letters();
    const t = await driven(b.getDocument(), { registry, size: { ring: [600, 500], lin: [360, 200] } });
    return { ...t, lin: t.ctx.viewports.get('lin') as Viewport, ring: t.ctx.layout.sectionIndex('ring') };
  }

  it("spans three nucleotides a residue, centered on the center residue's codon; a drag of 6 moves 2 residues", async () => {
    const t = await gene();
    const area = t.lin.live.current.width - (t.lin.axis?.insetX ?? 0);
    // Eight residues across from residue 2: the center residue's codon on the gene (a residue's fraction its codon's).
    t.lin.setView(area / 8, 2 * (area / 8), 'api');
    const c = centerOf(t.lin);
    const toGene = (r: number) => (3 * r < 34 ? 12 + 3 * r : 76 + 3 * r - 34);
    const g = new Float64Array(RING_WINDOW);
    expect(ringWindowFrame(t.ctx, t.ring, g)).toBe(true);
    const px = g[6];
    expect(((g[0] + g[1]) / 2 - g[2]) / px).toBeCloseTo(toGene(c.u), 0);
    expect((g[1] - g[0]) / px).toBeCloseTo((3 * area) / t.lin.live.current.cell, 0);
    const before = centerOf(t.lin).u;
    const l = t.ctx.layout.get('ring')?.current;
    const av = t.ctx.layout.axis('ring');
    if (!l || !av) throw new Error('no ring');
    const r = new Ring().from(l, av, 1);
    const p = new Float64Array(2);
    const q = new Float64Array(2);
    const mid = toGene(c.u);
    // On the rim (the ring is a ruler only).
    const rb = new Float64Array(2);
    ringBounds(l, rb);
    r.toSurface(av.blockX[0] + mid * px, rb[0], p);
    r.toSurface(av.blockX[0] + (mid + 6) * px, rb[0], q);
    expect(ringWindowAt(t.ctx, t.ring, { x: p[0], y: p[1] }, g)).toBe(WindowPart.pan);
    t.drag('ring', p[0], p[1], q[0], q[1]);
    expect(centerOf(t.lin).u - before).toBeCloseTo(2, 0);
  });
});
