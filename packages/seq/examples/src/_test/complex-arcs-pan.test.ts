/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The Complex's contacts as arcs while the view moves (core _spec/README.md "Render", The pairs layer, Arcs), the
// author's case: Examples › Complex › PAE and contacts › Contacts, the rail's window dragged across CDK2 (A), cyclin A
// (B) and the ATP chip. On a row of blocks an end past the window is taken at most `reach` of the window past the edge
// it lies beyond, so (1) a contact whose partner is off the window is an arc toward that edge, never a caret at the
// band's foot, never on the other side, its apex in the band; (2) panned a pixel at a time across the blocks' joints,
// every arc's points move a pixel or so, never jump (the scale counts the pairs within reach, so it does not snap); (3) a
// pan's frame equals a fresh renderer's at that window; (4) the contacts stay items at every window. The bugs: off a
// single strip every off-window end was a caret (the contacts went to ticks once the window held two blocks), and on
// one block the end was clamped in the section's units against the block's (a partner left of cyclin's block drawn
// rising off its right edge). Canvas and SVG draw one list: each arc's curve is in the SVG.

import { createSeqContext } from '@molstar/seqx-core/context/default';
import { DataStore, SeqRuntime } from '@molstar/seqx-core/data';
import {
  drawListsToSvg,
  type Frame,
  openModel,
  type PairsLay,
  PlacePoint,
  Renderer,
  renderConfig,
} from '@molstar/seqx-core/render';
import { describe, expect, it } from 'vitest';
import { nodeFixtureAssets } from '../node-assets.js';
import { ExampleViews } from '../views.js';

const Size = { width: 900, height: 400 };
/** The window's width in display units (CDK2 298, cyclin A's construct 260, the ATP chip 1). */
const Span = 150;
const Reach = 0.25;

async function opened() {
  const view = ExampleViews.find((v) => v.dataset === 'complex' && v.id === 'contacts')!;
  const ctx = createSeqContext({
    registry: [{ kind: 'format-runtime', runtime: SeqRuntime }, ...(view.registry?.({ delayMs: 0 }) ?? [])],
    ports: { assets: () => nodeFixtureAssets() },
  });
  const store = new DataStore();
  const o = await openModel(view.build(), store, {
    runtimes: ctx.runtimes,
    config: ctx.resolveConfig,
    assets: ctx.assets,
  });
  const renderer = () =>
    new Renderer({
      model: o.model,
      view: o.view!,
      data: store,
      tokens: ctx.tokens,
      text: ctx.text,
      types: ctx.types,
      layouts: ctx.layouts,
      sizings: ctx.sizings,
      config: renderConfig({}),
    });
  const r = renderer();
  let f = r.frame({ size: Size, settled: true });
  for (let i = 0; i < 6; i++) {
    store.want(f.wants);
    await store.whenSettled();
    f = r.frame({ size: Size, settled: true });
  }
  return { r, renderer };
}

const Opened = opened();

const at = (r: Renderer, start: number) =>
  r.frame({ size: Size, settled: true, viewport: { start, end: start + Span } });

function pairsLays(f: Frame): PairsLay[] {
  const out: PairsLay[] = [];
  for (const l of f.layout.sections[0].sectionLays) if (l.type.name === 'pairs') out.push(l.lay as PairsLay);
  return out;
}

const Pt = new PlacePoint();

/** Each end's display unit (the centre of its range on its system's block) of pair k: `[a, b]`. */
function unitsOf(f: Frame, lay: PairsLay, k: number): [number, number] {
  const sec = f.layout.sections[0];
  const d = lay.prep!.data!;
  const unit = (system: string, s: number, e: number) => {
    const b = sec.section.blocks!.find((b) => b.axis.system === system)!;
    return b.offset! + (s + e) / 2 - b.axis.ranges![0];
  };
  return [unit(lay.prep!.aSystem, d.aStart[k], d.aEnd[k]), unit(lay.prep!.bSystem, d.bStart[k], d.bEnd[k])];
}

/** Each arc's eight coordinates by its pair (a pair once on merged copies). */
function curves(f: Frame): Map<string, number[]> {
  const out = new Map<string, number[]>();
  for (const lay of pairsLays(f))
    for (let i = 0; i < lay.n; i++)
      out.set(`${lay.key}:${lay.pair[i]}`, Array.from(lay.curve.subarray(8 * i, 8 * i + 8)));
  return out;
}

/** Samples along a curve (a cubic's eight coordinates): x and y at t = i / 20. */
function samples(c: readonly number[]): number[] {
  const out: number[] = [];
  for (let i = 0; i <= 20; i++) {
    const t = i / 20;
    const u = 1 - t;
    for (let d = 0; d < 2; d++)
      out.push(u * u * u * c[d] + 3 * u * u * t * c[2 + d] + 3 * u * t * t * c[4 + d] + t * t * t * c[6 + d]);
  }
  return out;
}

/** The most any sample of a curve moved from another's. */
function moved(c: readonly number[], p: readonly number[]): number {
  const a = samples(c);
  const b = samples(p);
  return Math.max(...a.map((v, j) => Math.abs(v - b[j])));
}

describe('the Complex’s contacts as arcs while the view moves', () => {
  it('(1, 4) a partner off the window: an arc toward its edge, at most `reach` past it, in the band; never a caret', async () => {
    const { r } = await Opened;
    for (let start = 0; start + Span <= 560; start += 10) {
      const f = at(r, start);
      const sec = f.layout.sections[0];
      const ax = sec.axis!;
      const x0 = ax.x0;
      const x1 = ax.x0 + ax.width;
      const margin = Reach * ax.width;
      const s = sec.placement.strips;
      const svg = drawListsToSvg(f.lists);
      for (const lay of pairsLays(f)) {
        const where = `${lay.key} at ${start}`;
        // every end of these pairs is on a block of the section: none is a caret
        expect(lay.carets.length, where).toBe(0);
        expect(lay.densityN, where).toBe(0);
        let reachable = 0;
        for (let k = 0; k < lay.prep!.data!.n; k++) {
          const [ua, ub] = unitsOf(f, lay, k);
          const inA = ua >= ax.u0 && ua <= ax.u1;
          const inB = ub >= ax.u0 && ub <= ax.u1;
          if (inA || inB || ua < ax.u0 !== ub < ax.u0) reachable++;
        }
        // every pair with an end in view, or across it, is an arc
        expect(lay.n, where).toBe(reachable);
        for (let i = 0; i < lay.n; i++) {
          const c = lay.curve.subarray(8 * i, 8 * i + 8);
          const [ua, ub] = unitsOf(f, lay, lay.pair[i]);
          const lo = Math.min(ua, ub);
          const hi = Math.max(ua, ub);
          // the feet: the one in view where its unit is placed, the other past the edge its unit lies beyond
          const feet = [c[0], c[6]];
          for (const [u, x] of [
            [lo, feet[0]],
            [hi, feet[1]],
          ]) {
            if (u >= ax.u0 && u <= ax.u1 && sec.placement.place(u, Pt)) expect(x, where).toBeCloseTo(Pt.x, 6);
            else if (u < ax.u0) {
              expect(x, where).toBeLessThanOrEqual(x0 + 1e-6);
              expect(x, where).toBeGreaterThanOrEqual(x0 - margin - 1e-6);
            } else {
              expect(x, where).toBeGreaterThanOrEqual(x1 - 1e-6);
              expect(x, where).toBeLessThanOrEqual(x1 + margin + 1e-6);
            }
          }
          // the apex (the curve's middle) in the band and off the feet: never flat
          const apex = (c[1] + 3 * c[3] + 3 * c[5] + c[7]) / 8;
          let band = -1;
          for (let b = 0; b < s.bands; b++)
            if (c[1] >= s.bandY[b] + s.dy[0] && c[1] <= s.bandY[b] + s.bandH[b] + s.dy[0]) band = b;
          expect(band, where).toBeGreaterThanOrEqual(0);
          expect(apex, where).toBeGreaterThanOrEqual(s.bandY[band] + s.dy[0] - 1e-6);
          expect(c[1] - apex, where).toBeGreaterThanOrEqual(2);
        }
        // the SVG draws each arc's curve (canvas and SVG replay one list), to its 32-bit points
        const drawn = [...svg.matchAll(/d="M ([^"C]*) C ([^"C]*)"/g)].map((m) =>
          `${m[1]} ${m[2]}`.split(' ').map(Number),
        );
        for (let i = 0; i < lay.n; i++) {
          const c = lay.curve.subarray(8 * i, 8 * i + 8);
          expect(
            drawn.some((q) => q.every((v, j) => Math.abs(v - c[j]) < 0.02)),
            where,
          ).toBe(true);
        }
      }
    }
  });

  it('(2) panned a pixel at a time across the joints: every arc’s points move a pixel or so, never jump', async () => {
    const { r } = await Opened;
    const ax = at(r, 0).layout.sections[0].axis!;
    const step = Span / ax.width;
    const edges = [ax.x0, ax.x0 + ax.width];
    // across CDK2's end entering, then leaving, and the ATP chip entering (partners and the widest arcs coming and going)
    for (const [from, to] of [
      [140, 160],
      [290, 305],
      [400, 412],
    ]) {
      let prev: Map<string, number[]> | undefined;
      let prevUnits: number[] = [];
      for (let start = from; start <= to; start += step) {
        const f = at(r, start);
        const next = curves(f);
        // the residues' own move under the step (a chip entering takes its slot from the blocks): the arcs follow it
        const units: number[] = [];
        for (let u = Math.ceil(from); u < to + Span; u += 3)
          units.push(f.layout.sections[0].placement.place(u + 0.5, Pt) ? Pt.x : Number.NaN);
        let shift = 0;
        for (let j = 0; j < units.length; j++)
          if (Number.isFinite(units[j]) && Number.isFinite(prevUnits[j]))
            shift = Math.max(shift, Math.abs(units[j] - prevUnits[j]));
        const bound = Math.max(2, shift + 0.25);
        if (prev) {
          for (const [key, c] of next) {
            const p = prev.get(key);
            const where = `${key} at ${start.toFixed(2)}`;
            if (p) expect(moved(c, p), where).toBeLessThanOrEqual(bound);
            // a new arc enters at an edge: a foot there
            else
              expect(Math.min(...[c[0], c[6]].flatMap((x) => edges.map((e) => Math.abs(x - e)))), where).toBeLessThan(
                bound,
              );
          }
          for (const [key, p] of prev)
            if (!next.has(key))
              expect(Math.min(...[p[0], p[6]].flatMap((x) => edges.map((e) => Math.abs(x - e)))), key).toBeLessThan(
                bound,
              );
        }
        prev = next;
        prevUnits = units;
      }
    }
  });

  it('(3) a pan’s frame equals a fresh renderer’s full frame at that window', async () => {
    const { r, renderer } = await Opened;
    for (let start = 100; start <= 420; start += 37) {
      const panned = curves(at(r, start));
      const fresh = curves(at(renderer(), start));
      expect(panned, `at ${start}`).toEqual(fresh);
    }
  });
});
