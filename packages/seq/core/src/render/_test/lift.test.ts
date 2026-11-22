/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// A liftover (_spec/README.md "Render", A liftover; EC-212): the old assembly's genes drawn on the new one's blocks
// through a map read by region, its runs read by the new side (`to`).

import { SeqRuntime } from '@molstar/seqx-core/data';
import { loadGenome } from '@molstar/seqx-core/resolve';
import { beforeAll, describe, expect, it } from 'vitest';
import { loadBridges } from '../bridge-load.js';
import type { Frame } from '../renderer.js';
import type { RangesLay } from '../ranges.js';
import type { BlockPlacement } from '../layout/blocks.js';
import { settle } from './helpers.js';
import { Genes, imagesOf, liftSpec, NewChromosomes, NewSection, OldSection } from './lift-fixtures.js';

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
