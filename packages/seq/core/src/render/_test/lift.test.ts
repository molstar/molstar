/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// A liftover (_spec/README.md "Render", A liftover; EC-212): the old assembly's genes drawn on the new one's blocks
// through a map read by region, its runs read by the new side (`to`).

import { SeqRuntime } from '@molstar/seqx-core/data';
import { loadGenome } from '@molstar/seqx-core/resolve';
import { beforeAll, describe, expect, it, vi } from 'vitest';
import { TracksLayout } from '../index.js';
import { RibbonsType } from '../pairs-ribbons.js';
import { loadBridges } from '../bridge-load.js';
import type { Frame } from '../renderer.js';
import type { RangesLay } from '../ranges.js';
import type { BlockPlacement } from '../layout/blocks.js';
import { settle } from './helpers.js';
import { Genes, imagesOf, liftSpec, NewChromosomes, NewSection, OldSection, Runs } from './lift-fixtures.js';
import { regionProjector } from '../lift-project.js';

beforeAll(async () => {
  await loadGenome();
  await loadBridges();
});

const Size = { width: 900, height: 400 };
const chrOf = (system: string) => /:name (chr\d)\)/.exec(system)?.[1] ?? system;

/** The rows a section's `genes` layer drew, by block: key, start, end, sorted. */
function drawn(f: Frame, section: string): Record<string, [string, number, number][]> {
  const sec = f.layout.sections.find((s) => s.section.address === section)!;
  const strips = (sec.placement as BlockPlacement).strips;
  const bands = sec.lays.length / strips.n;
  const out: Record<string, [string, number, number][]> = {};
  for (let i = 0; i < strips.n; i++) {
    const rows: [string, number, number][] = [];
    for (let k = 0; k < bands; k++)
      for (const l of sec.lays[i * bands + k] ?? []) {
        if (!l.key.endsWith('/genes')) continue;
        const r = (l.lay as RangesLay).rows;
        for (let j = 0; j < r.n; j++) rows.push([String(r.keys[j]), r.starts[j], r.ends[j]]);
      }
    out[chrOf(strips.axes[i].system)] = rows.sort((x, y) => x[1] - y[1] || x[0].localeCompare(y[0]));
  }
  return out;
}

/** Every gene's images through the runs, by new chromosome: the truth. */
function truth(): Record<string, [string, number, number][]> {
  const out: Record<string, [string, number, number][]> = {};
  for (const c of NewChromosomes) out[c.name] = [];
  Genes.id.forEach((id, i) => {
    for (const m of imagesOf(Genes.target[i], Genes.start[i], Genes.end[i])) out[m.b].push([id, m.start, m.end]);
  });
  for (const k of Object.keys(out)) out[k].sort((x, y) => x[1] - y[1] || x[0].localeCompare(y[0]));
  return out;
}

describe('a layer lifted through the map', () => {
  it('draws each item’s images through the runs: many to one, one to many, a deletion, an inversion (EC-20, EC-21)', async () => {
    const s = await settle(liftSpec(), [SeqRuntime]);
    const f = await s.frame({ size: Size });
    expect(drawn(f, NewSection)).toEqual(truth());
    // the old section draws the genes as they are
    expect(Object.values(drawn(f, OldSection)).flat().length).toBe(Genes.id.length);
  });

  it('a part that maps nowhere is not drawn; an item wholly deleted is not drawn at all (EC-21)', async () => {
    const s = await settle(liftSpec(), [SeqRuntime]);
    const f = await s.frame({ size: Size });
    const rows = Object.values(drawn(f, NewSection)).flat();
    expect(rows.filter((r) => r[0] === 'lost')).toEqual([]);
    expect(rows.filter((r) => r[0] === 'half')).toEqual([['half', 20_000, 25_000]]);
    expect(rows.filter((r) => r[0] === 'twice').map((r) => r.slice(1))).toEqual([
      [60_000 + 2_000, 60_000 + 8_000],
      [60_000 + 2_000, 60_000 + 8_000],
    ]);
  });

  it('reads the runs by the new side within the window: zoomed in, only the images there (EC-212)', async () => {
    const s = await settle(liftSpec(), [SeqRuntime]);
    // the new section shows chr1 then chr2 end to end: a window over chr1's 45–70 kb
    const f = await s.frame({ size: Size, viewports: { [NewSection]: { start: 45_000, end: 70_000 } } });
    const rows = drawn(f, NewSection).chr1 ?? [];
    expect(rows).toEqual(truth().chr1.filter((r) => r[1] < 70_000 && r[2] > 45_000));
    // the reads by side `b` alone on the new section's blocks
    const wants = f.wants.filter((w) => w.source === '/data/chain' && w.layer?.endsWith('/new/genes/genes'));
    expect(wants.length).toBeGreaterThan(0);
    expect(wants.every((w) => w.side === 'b' && w.items)).toBe(true);
  });
});

describe('a pan of the new section', () => {
  it('lays out that section and the bridge only; a frame with nothing moved, neither (EC-212)', async () => {
    const s = await settle(liftSpec(), [SeqRuntime]);
    const at = (k: number) => ({
      size: Size,
      viewports: { [NewSection]: { start: 1_000 + k, end: 120_000 + k } },
    });
    await s.frame(at(0));
    await s.frame(at(0));
    const place = vi.spyOn(TracksLayout, 'place');
    const ribbons = vi.spyOn(RibbonsType, 'layout');
    s.renderer.frame(at(5));
    expect(place.mock.calls.map((c) => c[0].address)).toEqual([NewSection]);
    expect(ribbons).toHaveBeenCalledTimes(1);
    place.mockClear();
    ribbons.mockClear();
    s.renderer.frame({ ...at(5), emphasis: { hover: { layer: '', section: '', system: '' } } });
    expect(place).not.toHaveBeenCalled();
    expect(ribbons).not.toHaveBeenCalled();
    place.mockRestore();
    ribbons.mockRestore();
  });
});

/** Each lifted row of the new section: its key, extent, thick part and strand as carried; each item's chevrons. */
function carried(f: Frame) {
  const sec = f.layout.sections.find((s) => s.section.address === NewSection)!;
  const strips = (sec.placement as BlockPlacement).strips;
  const bands = sec.lays.length / strips.n;
  const out: { chr: string; key: string; start: number; end: number; ts: number; te: number; strand: number }[] = [];
  const chevrons: Record<string, number> = {};
  for (let i = 0; i < strips.n; i++)
    for (let k = 0; k < bands; k++)
      for (const l of sec.lays[i * bands + k] ?? []) {
        if (!l.key.endsWith('/genes')) continue;
        const lay = l.lay as RangesLay;
        const r = lay.rows;
        for (let j = 0; j < r.n; j++) {
          const col = (name: string) => r.fields[j]?.find((c) => c.name === name)?.column as ArrayLike<number>;
          out.push({
            chr: chrOf(strips.axes[i].system),
            key: String(r.keys[j]),
            start: r.starts[j],
            end: r.ends[j],
            ts: col('thick_start')[r.row[j]],
            te: col('thick_end')[r.row[j]],
            strand: col('strand')[r.row[j]],
          });
        }
        for (let c = 0; c < lay.cn; c++) chevrons[String(r.keys[lay.citem[c]])] = lay.cdir[c];
      }
  return { rows: out.sort((a, b) => a.chr.localeCompare(b.chr) || a.start - b.start), chevrons };
}

describe('a lifted gene’s own coordinates', () => {
  it('its thick part carried and clipped to each image, its strand turned by an inverted run (EC-21, EC-212)', async () => {
    const s = await settle(liftSpec(), [SeqRuntime]);
    const f = await s.frame({ size: Size });
    const { rows: got, chevrons } = carried(f);
    const want: typeof got = [];
    Genes.id.forEach((key, i) => {
      const thick = imagesOf(Genes.target[i], Genes.thick_start[i], Genes.thick_end[i]);
      for (const m of imagesOf(Genes.target[i], Genes.start[i], Genes.end[i])) {
        let lo = Number.POSITIVE_INFINITY;
        let hi = Number.NEGATIVE_INFINITY;
        for (const t of thick)
          if (t.b === m.b && Math.min(t.end, m.end) > Math.max(t.start, m.start)) {
            lo = Math.min(lo, Math.max(t.start, m.start));
            hi = Math.max(hi, Math.min(t.end, m.end));
          }
        const strand = (Genes.strand[i] === '-' ? -1 : 1) * m.dir;
        const thin = !(hi > lo);
        want.push({
          chr: m.b,
          key,
          start: m.start,
          end: m.end,
          ts: thin ? m.start : lo,
          te: thin ? m.start : hi,
          strand,
        });
      }
    });
    expect(got).toEqual(want.sort((a, b) => a.chr.localeCompare(b.chr) || a.start - b.start));
    // the inverted gene: its thick part 56–59 kb on old chr1 lands at 51–54 kb, read the other way
    expect(got.find((r) => r.key === 'inverted')).toMatchObject({ ts: 51_000, te: 54_000, strand: -1 });
    // the half-deleted gene: its thick part only where it lands
    expect(got.find((r) => r.key === 'half')).toMatchObject({ start: 20_000, end: 25_000, ts: 20_000, te: 24_000 });
    // the chevrons by the strand as carried: the inverted `+` gene points left, a straight `+` right, a `-` left
    expect(chevrons.inverted).toBe(2);
    expect(chevrons.half).toBe(1);
    expect(chevrons.straight).toBe(2);
  });

  it('rows whose ends differ in length are no runs, said (`seq.render.lift.run`)', async () => {
    const unequal = {
      a: 'chr2',
      from: 40_000,
      b: 'chr2',
      to: 70_000,
      length: 5_000,
      bLength: 4_000,
      strand: 1 as const,
    };
    const s = await settle(liftSpec([...Runs, unequal]), [SeqRuntime]);
    const f = await s.frame({ size: Size });
    const said = f.layout.issues.filter((i) => i.code === 'seq.render.lift.run');
    expect(said.length).toBeGreaterThan(0);
    expect(said[0].message).toContain('1 rows');
  });
});

describe('places, through the map (EC-20)', () => {
  it('a duplicated gene in two places from either side; one split by an insertion longer than it in one', async () => {
    const s = await settle(liftSpec(), [SeqRuntime]);
    await s.frame({ size: Size });
    const p = regionProjector(s.store);
    const old = (c: string) => `/data/old/(sequence :name ${c})`;
    const neu = (c: string) => `/data/new/(sequence :name ${c})`;
    const two = '2 places: chr1 62,001–68,000, chr2 62,001–68,000';
    expect(p.places!(s.model, s.view, old('chr1'), 72_000, 78_000)).toBe(two);
    // hovered on one image: the item's places, itself among them
    expect(p.places!(s.model, s.view, neu('chr2'), 62_000, 68_000)).toBe(two);
    // 25–30 kb and 45–50 kb of new chr2: one colinear chain, one place
    expect(p.places!(s.model, s.view, old('chr2'), 15_000, 25_000)).toBeUndefined();
    // the other image lit from one, through the old gene
    const out: number[] = [];
    expect(p.project(s.model, neu('chr2'), [62_000, 68_000], neu('chr1'), out)).toBe(true);
    expect(out).toEqual([62_000, 68_000]);
  });
});
