/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// A genome's row of blocks (_spec/README.md "Render", Blocks): a pan across a block gap places
// the row (its gaps, runs and groups, the window's inverse) with nothing allocated once warm; a warm draw of the row
// with its runs and a group's bracket, and a warm hover on it, allocate nothing; a pan frame across the gap (a layout
// and a draw) allocates a bounded size. By hand, never in CI: `pnpm -C packages/seq/core test:perf`.

import { quickCase } from '@molstar/test-harness/quick';
import { beforeAll, describe, expect, inject, it } from 'vitest';
import { loadGenome } from '@molstar/seqx-core/resolve';
import { type BlockPlacement, drawView, type FrameInput, type Hit, hitView } from '../index.js';
import { categoryColumn } from '@molstar/core/repr/tree';
import { type Coverage, type PairsList, type SeqData, SeqRuntime, seqToValue } from '@molstar/seqx-core/data';
import { assembly, createSeqSpec, genome } from '@molstar/seqx-core/lang';
import { genomeSpec, toyGenome } from './genome-fixtures.js';
import { settle } from './helpers.js';
import { reusedHover } from './scenarios.js';

const Size = { width: 1200, height: 300 };
const Steps = 64;

beforeAll(async () => {
  await loadGenome();
});

/** A toy genome of 25 chromosomes and 400 contigs, a window across chr3's end into chr4 (the gap in its middle). */
async function genomeRow() {
  const s = await settle(genomeSpec(toyGenome(25, 400), { reversed: ['chr4'], scaled: { chr5: 2 } }));
  const fit = await s.frame({ size: Size, dpr: 2 });
  const sec = fit.layout.sections[0].section;
  const gap = sec.blocks![3].offset!;
  const w = 6e6;
  const inputs: FrameInput[] = [];
  for (let k = 0; k < Steps; k++) {
    const start = gap - w / 2 + ((k - Steps / 2) * w) / (2 * Steps);
    inputs.push({ size: Size, dpr: 2, viewport: { start, end: start + w } });
  }
  for (const i of inputs.slice(0, 4)) await s.frame(i);
  // the fit again, last: a frame's placement is kept and reused by the next
  return { s, fit: await s.frame({ size: Size, dpr: 2 }), inputs };
}

describe('a genome’s row of blocks: pans across a gap, draws and hovers', () => {
  it('the row placed for a pan across a block gap: 0 bytes warm (B6)', async ({ task }) => {
    const { s, inputs } = await genomeRow();
    const f = s.renderer.frame(inputs[0]);
    const p = f.layout.sections[0].placement as BlockPlacement;
    const row = p.row;
    const windows = inputs.map((i) => i.viewport!);
    const { r } = quickCase(task, 'genome row pan across a gap', (i) => {
      const v = windows[i % Steps];
      row.lay(v.start, v.end);
      row.snap(2);
      // the pointer's unit back from its px: the inverse a pan reads
      row.unitAt(row.x0 + row.width / 2);
    });
    console.info(`genome row pan across a gap: ${r.bytesPerStep.toFixed(0)} bytes a step (${r.note})`);
    expect(r.bytesPerStep, r.note).toBeLessThanOrEqual(inject('perfThreshold'));
  });

  it('a pan frame across the gap (layout and draw): a bounded size', async ({ task }) => {
    const { s, inputs } = await genomeRow();
    const { r } = quickCase(task, 'genome row pan frame', (i) => void s.renderer.frame(inputs[i % Steps]), {
      warmSteps: 300,
      steps: 120,
      windows: 3,
      maxWarmSteps: 2000,
    });
    console.info(`genome row pan frame: ${r.bytesPerStep.toFixed(0)} bytes a frame (${r.note})`);
    // measured 31.7 KB (2026-11-21): the tracks' lays of two blocks, the ruler's ticks and labels, the key; 31.6 KB
    // with the block labels' strings made once per resolved block (two labels in this window: about 0.1 KB of it)
    expect(r.bytesPerStep, r.note).toBeLessThan(48_000);
  });

  it('the whole genome at the fit, its runs and labels: draw', async ({ task }) => {
    const { fit } = await genomeRow();
    expect((fit.layout.sections[0].placement as BlockPlacement).runN).toBeGreaterThan(0);
    const { r } = quickCase(task, 'genome row draw', () => drawView(fit.layout, fit.lists), {});
    console.info(`genome row draw: ${r.bytesPerStep.toFixed(0)} bytes a frame (${r.note})`);
    expect(r.bytesPerStep, r.note).toBeLessThanOrEqual(inject('perfThreshold'));
  });

  it('the whole genome at the fit: hover', async ({ task }) => {
    const { s, fit } = await genomeRow();
    const p = fit.layout.sections[0].placement as BlockPlacement;
    const hits: Hit[] = [];
    for (let i = 0; i < p.strips.n; i++) {
      const a = p.strips.axes[i];
      const h = hitView(fit.layout, a.x0 + a.width / 2, p.strips.bandY[i * p.strips.bands] + 4);
      if (h) hits.push(h);
    }
    const hover = reusedHover(hits);
    const at: FrameInput = { size: Size, dpr: 2, emphasis: hover.emphasis };
    const { r } = quickCase(task, 'genome row hover', (i) => {
      hover.set(i);
      s.renderer.frame(at);
    });
    console.info(`genome row hover: ${r.bytesPerStep.toFixed(0)} bytes a frame (${r.note})`);
    expect(r.bytesPerStep, r.note).toBeLessThanOrEqual(inject('perfThreshold'));
  });
});

/** 25 chromosomes, a features layer over the assembly (rows split by `target`) and links over it as a ring's chords. */
async function genomeData() {
  const chromosomes = toyGenome(25, 0);
  const b = createSeqSpec({ title: 'data over 25 blocks' });
  const g = b.data.assembly('g', assembly('toy', chromosomes, { order: 'source' } as never));
  const start: number[] = [];
  const end: number[] = [];
  const target: string[] = [];
  chromosomes.forEach((c, i) => {
    for (let k = 0; k < 40; k++) {
      start.push(Math.floor((k * c.length) / 40));
      end.push(Math.floor((k * c.length) / 40) + 400_000);
      target.push(i % 2 ? c.name.slice(3) : c.name);
    }
  });
  const features = b.data.ranges('features', { target: g, rows: { start, end, target } } as never);
  // 2,000 links between the chromosomes, each listed on both its ends' chromosomes, as a `seq` source over them
  const n = 2000;
  const ends = Array.from({ length: n }, (_, i) => {
    const a = i % chromosomes.length;
    const b2 = (i * 7) % chromosomes.length;
    return {
      a,
      b: b2,
      as: (i * 7919) % (chromosomes[a].length - 100),
      bs: (i * 104729) % (chromosomes[b2].length - 100),
    };
  });
  const pairs: SeqData['pairs'] = chromosomes.map((c, k) => {
    const mine = ends.filter((e) => e.a === k || e.b === k);
    const list: PairsList = {
      form: 'list',
      aSystem: categoryColumn(mine.map((e) => chromosomes[e.a].name)),
      bSystem: categoryColumn(mine.map((e) => chromosomes[e.b].name)),
      aStart: Int32Array.from(mine, (e) => e.as),
      aEnd: Int32Array.from(mine, (e) => e.as + 50),
      bStart: Int32Array.from(mine, (e) => e.bs),
      bEnd: Int32Array.from(mine, (e) => e.bs + 50),
    };
    return {
      request: { system: c.name, start: 0, end: c.length, level: 1, shapes: ['pairs'] },
      answer: { ok: true, tables: { pairs: list } },
    };
  });
  const coverage: Coverage = {
    systems: chromosomes.map((c) => ({ name: c.name, length: c.length })),
    shapes: chromosomes.map((c) => ({ system: c.name, shape: 'pairs' as const })),
    levels: [1],
  };
  const src = b.data.source('links-src', { format: 'seq', value: seqToValue({ coverage, pairs }) as never });
  const links = b.data.pairs('links', { a: g, source: src } as never);
  b.views.view('v', (v) =>
    genome(v, 'genome', g, {}, (sec) => {
      sec.track('features').layer('features', { type: 'ranges', data: features });
      sec.track('links').layer('links', { type: 'pairs', data: links });
    }),
  );
  const s = await settle(b.build(), [SeqRuntime]);
  const input: FrameInput = { size: Size, dpr: 2 };
  for (let i = 0; i < 4; i++) await s.frame(input);
  return { s, input };
}

describe('data over an assembly: 25 blocks', () => {
  it('a warm frame of a features and a pairs layer over 25 blocks (read by sequence): 0 bytes', async ({ task }) => {
    const { s, input } = await genomeData();
    const f = s.renderer.frame(input);
    expect((f.layout.sections[0].placement as BlockPlacement).strips.n).toBe(25);
    expect((f.layout.sections[0].sectionLays[0].lay as unknown as { n: number }).n).toBeGreaterThan(100);
    const { r } = quickCase(task, 'genome data warm frame', () => void s.renderer.frame(input));
    console.info(`genome data warm frame: ${r.bytesPerStep.toFixed(0)} bytes a frame (${r.note})`);
    expect(r.bytesPerStep, r.note).toBeLessThanOrEqual(inject('perfThreshold'));
  });

  it('its draw: 0 bytes', async ({ task }) => {
    const { s, input } = await genomeData();
    const f = s.renderer.frame(input);
    const { r } = quickCase(task, 'genome data draw', () => drawView(f.layout, f.lists), {});
    console.info(`genome data draw: ${r.bytesPerStep.toFixed(0)} bytes a frame (${r.note})`);
    expect(r.bytesPerStep, r.note).toBeLessThanOrEqual(inject('perfThreshold'));
  });
});
