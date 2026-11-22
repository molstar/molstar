/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { DarkTokens, LightTokens, type PairsLay, type SectionLayout } from '@molstar/seqx-core/render';
import { loadMembers } from '@molstar/seqx-core/render/members-hooks';
import { loadGenome } from '@molstar/seqx-core/resolve';
import { blocksOf, Inversion, ItemsPitch, PeachToGrape, Translocation } from '@molstar/seqx-examples/across/genomes';
import {
  ContactBins,
  contactValues,
  DotsSection,
  Hotspots,
  kinshipValues,
  Samples,
  SamplesSystem,
} from '@molstar/seqx-examples/across/grids';
import { AcrossViews } from '@molstar/seqx-examples/across/views';
import { describerFor } from '@molstar/seqx-viewer/describe';
import type { SeqViewer } from '@molstar/seqx-viewer/viewer';
import { afterEach, beforeAll, describe, expect, it, vi } from 'vitest';
import { openExample } from '../../../examples/_test/draw-checks.js';

// Features/Across sections, the grids, headless on a fake clock as the stories compose them (a viewer of the view, its
// SVG export light and dark): peach against grape as a dot plot over two sections of blocks (EC-209: the inversion an
// anti-diagonal, the translocation in another block pair's cell), contacts between two chromosomes' sections (EC-223),
// and a kinship grid on categorical axes beside its genotype matrix (EC-35).

let viewer: SeqViewer | undefined;

// the genome and categorical chunks landed first: their imports are real I/O, never on the fake clock
beforeAll(async () => {
  await loadGenome();
  await loadMembers();
});

afterEach(() => {
  viewer?.retire();
  viewer = undefined;
  vi.useRealTimers();
});

const view = (id: string) => AcrossViews.find((v) => v.id === id)!;

async function open(id: string) {
  viewer = await openExample(view(id), { delayMs: 300, wait: 30_000 });
  viewer.draw();
  return viewer;
}

/** Both themes' exports settled, with no issue, and different. */
async function bothThemes(v: SeqViewer) {
  const out: string[] = [];
  for (const dark of [false, true]) {
    const r = v.renderSvg({ idPrefix: dark ? 'gd-' : 'gl-', tokens: dark ? DarkTokens : LightTokens });
    await vi.advanceTimersByTimeAsync(4000);
    const svg = await r;
    expect(svg.settled).toBe(true);
    expect(svg.issues.map((i) => `${i.code}: ${i.message}`)).toEqual([]);
    out.push(svg.svg);
  }
  expect(out[0]).not.toBe(out[1]);
  return out[0];
}

const sectionOf = (v: SeqViewer, address: string): SectionLayout =>
  v.frame!.layout.sections.find((s) => s.section.address === address)!;
const named = (v: SeqViewer, system: string) => v.view.model.value.systems.get(system)?.name;
const card = (v: SeqViewer) => describerFor({ ...v, layout: () => v.frame?.layout }).describe(v.hover.current)!;

describe('Features/Across sections, the grids, headless', () => {
  it('peach against grape as a dot plot (EC-209): a block each, the inversion an anti-diagonal, the translocation in Pp03 × chr4', async () => {
    const v = await open('across-dot-plot');
    expect(v.layoutIssues.value).toEqual([]);
    const s = sectionOf(v, DotsSection);
    const strips = s.placement.strips;
    expect(strips.n).toBe(4 * 5);
    const blocks = blocksOf(PeachToGrape, ItemsPitch, 7);
    let drawn = 0;
    for (let i = 0; i < strips.n; i++) {
      const lay = s.lays[i][0]?.lay as PairsLay | undefined;
      const a = named(v, strips.axes[i].system);
      const b = named(v, strips.ys[i]!.system);
      const mine = blocks.a.filter((x, k) => x === a && blocks.b[k] === b).length;
      expect((lay?.nd ?? 0) + (lay?.ns ?? 0), `${a} × ${b}`).toBe(mine);
      drawn += mine;
      if (!lay) continue;
      for (let k = 0; k < lay.ns; k++) {
        const [x0, y0, x1, y1] = lay.segs.subarray(5 * k, 5 * k + 4);
        const against = y1 < y0 === x1 > x0;
        // the inversion's blocks alone run against the diagonal (y falls as x grows)
        const d = lay.prep!.data!;
        const p = lay.segPair[k];
        const inverted = a === 'Pp02' && d.aStart[p] >= Inversion[0] && d.aEnd[p] <= Inversion[1];
        expect(against, `${a}:${d.aStart[p]}`).toBe(inverted);
      }
    }
    expect(drawn).toBe(blocks.a.length);
    // the translocation: Pp03's tail on grape's chr4
    expect(blocks.a.some((x, k) => x === 'Pp03' && blocks.b[k] === 'chr4' && blocks.aStart[k] >= Translocation)).toBe(
      true,
    );
    await bothThemes(v);
  });

  it('contacts between two chromosomes’ sections (EC-223): each hotspot under its cell, the strongest values', async () => {
    const v = await open('across-contacts');
    expect(v.layoutIssues.value).toEqual([]);
    const s = sectionOf(v, '/views/main/contacts');
    const strips = s.placement.strips;
    expect(strips.n).toBe(1);
    const values = contactValues();
    for (const [i, j] of Hotspots) {
      strips.axes[0].cell(i);
      const x = (strips.axes[0].sx0 + strips.axes[0].sx1) / 2;
      strips.ys[0]!.cell(j);
      const y = (strips.ys[0]!.sx0 + strips.ys[0]!.sx1) / 2;
      const h = v.hit(x, y)!;
      expect([h.position, h.yPosition, h.block, h.yBlock]).toEqual([i, j, 0, 0]);
      expect(values[i * ContactBins.b + j]).toBeGreaterThan(0.9);
    }
    await bothThemes(v);
  });

  it('kinship beside the genotypes (EC-35): members by label, families on the diagonal, a label a sample', async () => {
    const v = await open('across-kinship');
    expect(v.layoutIssues.value).toEqual([]);
    const s = sectionOf(v, '/views/main/kinship');
    const a = s.placement.strips.axes[0];
    const ay = s.placement.strips.ys[0]!;
    const sorted = [...Samples].sort();
    expect(a.members!.order.length).toBe(Samples.length);
    expect([...a.members!.order].map((m) => Samples[m])).toEqual(sorted);
    expect([...ay.members!.order].map((m) => Samples[m])).toEqual(sorted);
    // the truth the story shows: a family's members more alike than strangers
    const k = kinshipValues();
    const n = Samples.length;
    const fam = (i: number) => Samples[i].slice(0, 2);
    let within = 0;
    let wn = 0;
    let between = 0;
    let bn = 0;
    for (let i = 0; i < n; i++)
      for (let j = 0; j < n; j++)
        if (i === j) continue;
        else if (fam(i) === fam(j)) (within += k[i * n + j]), wn++;
        else (between += k[i * n + j]), bn++;
    expect(within / wn).toBeGreaterThan(between / bn + 0.2);
    // a cell's card names both members
    const cx = (a.xOfUnit(0) + a.xOfUnit(1)) / 2;
    const cy = (ay.xOfUnit(1) + ay.xOfUnit(2)) / 2;
    // the hover model names a section by its address
    v.hover.set({ ...v.hit(cx, cy)!, section: s.section.address }, 'test');
    expect(card(v).title).toBe(`${sorted[0]} × ${sorted[1]}`);
    // a label is its sample, a row of the samples' system
    const label = v.hit(a.x0 - 20, cy)!;
    expect([label.target, label.system, Samples[label.position!]]).toEqual(['row', SamplesSystem, sorted[1]]);
    const svg = await bothThemes(v);
    for (const name of Samples) expect(svg).toContain(name);
  });
});
