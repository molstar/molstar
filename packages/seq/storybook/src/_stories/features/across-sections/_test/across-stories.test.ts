/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { DarkTokens, LightTokens } from '@molstar/seqx-core/render';
import { loadBridges } from '@molstar/seqx-core/render/bridge-load';
import { loadGenome } from '@molstar/seqx-core/resolve';
import {
  AcrossViews,
  GrapeSection,
  PeachBelowSection,
  PeachSection,
  PlumSection,
} from '@molstar/seqx-examples/across/views';
import {
  BinsPitch,
  type Blocks,
  blocksOf,
  Gap,
  Inversion,
  ItemsPitch,
  PeachToGrape,
  Translocation,
} from '@molstar/seqx-examples/across/genomes';
import { describerFor } from '@molstar/seqx-viewer/describe';
import type { SeqViewer } from '@molstar/seqx-viewer/viewer';
import { afterEach, beforeAll, describe, expect, it, vi } from 'vitest';
import { openExample } from '../../../examples/_test/draw-checks.js';
import { openZoomed, openZoomedGrape } from '../openings.js';

// Features/Across sections, headless on a fake clock as the stories compose them (a viewer of the view, its SVG export
// light and dark): two genomes at the fit (an inversion crossed, a translocation, a gap; EC-206, EC-218), zoomed out to
// bins (EC-129), one side zoomed in, a stub, the tiles arriving a side at a time (the pending hatch, EC-209), a reader
// without bridges, one genome against itself with `mirror`, three genomes chained (EC-235).

let viewer: SeqViewer | undefined;

// the genome and bridge chunks landed first: their imports are real I/O, never on the fake clock
beforeAll(async () => {
  await loadGenome();
  await loadBridges();
  await import('@molstar/seqx-viewer/bridge-view');
});

afterEach(() => {
  viewer?.retire();
  viewer = undefined;
  vi.useRealTimers();
});

/** The read's delay (`a` side), ms: the `b` side's is twice it. */
const DelayMs = 300;

const view = (id: string) => AcrossViews.find((v) => v.id === id)!;

async function open(id: string, wait = 30_000) {
  viewer = await openExample(view(id), { delayMs: DelayMs, wait });
  viewer.draw();
  return viewer;
}

async function exported(v: SeqViewer, dark: boolean) {
  const r = v.renderSvg({ idPrefix: dark ? 'ad-' : 'al-', tokens: dark ? DarkTokens : LightTokens });
  await vi.advanceTimersByTimeAsync(4000);
  return r;
}

const bridges = (v: SeqViewer) => v.frame!.layout.bridges?.list ?? [];

/** Both themes' exports settled, with no issue, a bridge group each, and different. */
async function bothThemes(v: SeqViewer, groups = 1) {
  const out: string[] = [];
  for (const dark of [false, true]) {
    const svg = await exported(v, dark);
    expect(svg.settled).toBe(true);
    expect(svg.issues.map((i) => `${i.code}: ${i.message}`)).toEqual([]);
    expect(svg.svg.match(/data-bridge=/g)?.length ?? 0).toBe(groups);
    out.push(svg.svg);
  }
  expect(out[0]).not.toBe(out[1]);
}

const named = (v: SeqViewer, system: string) => v.view.model.value.systems.get(system)?.name;

/** Whether a section's window shows a position of a chromosome (a block's place over its window). */
function shows(v: SeqViewer, section: string, chr: string, p: number): boolean {
  const s = v.frame!.layout.sections.find((x) => x.section.address === section)!;
  const strips = s.placement.strips;
  for (let i = 0; i < strips.n; i++) {
    const axis = strips.axes[i];
    if (named(v, axis.system) !== chr) continue;
    const vis = axis.visible([]);
    for (let j = 0; j + 1 < vis.length; j += 2) if (p >= vis[j] && p < vis[j + 1]) return true;
  }
  return false;
}

/**
 * Every block anchored in either side's window drawn exactly once (its item, or a bin pair of its strand holding both
 * its starts), none twice (EC-209).
 */
function expectDrawnOnce(v: SeqViewer, bl: Blocks) {
  const ribbons = bridges(v)[0].ribbons();
  const count = new Array<number>(bl.a.length).fill(0);
  for (const r of ribbons) {
    const a = named(v, r.aSystem);
    const b = named(v, r.bSystem);
    for (let j = 0; j < bl.a.length; j++) {
      if (bl.a[j] !== a || bl.b[j] !== b || r.reversed !== bl.strand[j] < 0) continue;
      const as = bl.aStart[j];
      const bs = bl.bStart[j];
      if (
        r.bin === undefined
          ? r.aStart === as && r.bStart === bs
          : as >= r.aStart && as < r.aEnd && bs >= r.bStart && bs < r.bEnd
      )
        count[j]++;
    }
  }
  const wrong: string[] = [];
  for (let j = 0; j < bl.a.length; j++) {
    const shown = shows(v, PeachSection, bl.a[j], bl.aStart[j]) || shows(v, GrapeSection, bl.b[j], bl.bStart[j]);
    if (count[j] > 1 || (shown && count[j] !== 1))
      wrong.push(`${bl.a[j]}:${bl.aStart[j]} ${bl.b[j]}:${bl.bStart[j]} ×${count[j]}`);
  }
  expect(wrong.slice(0, 5)).toEqual([]);
}

describe('Features/Across sections, headless', () => {
  it('two genomes at the fit (EC-206, EC-218): a ribbon a block, the inversion crossed, the translocation to chr4, no ribbon over the gap', async () => {
    const v = await open('across-fit');
    const [b] = bridges(v);
    const ribbons = b.ribbons();
    expect(ribbons.length).toBeGreaterThan(200);
    const on = (r: (typeof ribbons)[number]) => [named(v, r.aSystem), named(v, r.bSystem)];
    // the inversion: Pp02's blocks inside it crossed, the others not
    const chr2 = ribbons.filter((r) => on(r)[0] === 'Pp02');
    expect(chr2.filter((r) => r.reversed).every((r) => r.aStart >= Inversion[0] && r.aEnd <= Inversion[1])).toBe(true);
    expect(chr2.some((r) => r.reversed)).toBe(true);
    // the translocation: Pp03's tail on grape's chr4
    expect(ribbons.some((r) => on(r)[0] === 'Pp03' && on(r)[1] === 'chr4' && r.aStart >= Translocation)).toBe(true);
    // the gap: nothing over Pp01's lost stretch
    expect(ribbons.some((r) => on(r)[0] === 'Pp01' && r.aEnd > Gap[0] && r.aStart < Gap[1])).toBe(false);
    await bothThemes(v);
  });

  it('zoomed out (EC-129): past the gate, bin pairs, never the items', async () => {
    const v = await open('across-bins');
    const [b] = bridges(v);
    const lay = b.lays[0].lay as unknown as { reads: { level?: number }[]; n: number };
    expect(lay.n).toBeGreaterThan(0);
    expect(lay.reads.some((r) => (r.level ?? 1) > 1)).toBe(true);
    // a bin's card: its count and its bins' width, never an aligned length
    const r = b.ribbons().find((x) => x.bin !== undefined && x.ax1 - x.ax0 > 1 && x.ax0 > b.box.x)!;
    v.hover.set(v.hit((r.ax0 + r.ax1) / 2, b.box.y + 1), 'test');
    const card = describerFor({ ...v, layout: () => v.frame?.layout }).describe(v.hover.current)!;
    expect(card.text).toMatch(/pairs? · .* bins · zoom in to see them/);
    expect(card.rows.map((row) => row.label)).not.toContain('Aligned');
    await bothThemes(v);
  });

  it('one side zoomed in: peach in Pp02’s 300 kb, its ribbons fanning to grape at the fit', async () => {
    const v = await open('across-zoomed');
    await openZoomed(v);
    await vi.advanceTimersByTimeAsync(4000);
    v.draw();
    const w = v.view.get(PeachSection)!;
    const g = v.view.get(GrapeSection)!;
    expect(w.end - w.start).toBeLessThan(400_000);
    expect(g.start).toBe(0);
    const lay = bridges(v)[0].lays[0].lay as unknown as { reads: unknown[]; readSide: number[] };
    // each side its own read: peach's window by peach, the rest of grape's blocks by grape
    expect(lay.readSide).toEqual([0, 1]);
    const ribbons = bridges(v)[0].ribbons();
    expect(ribbons.filter((r) => named(v, r.aSystem) === 'Pp02').length).toBeGreaterThan(5);
    expectDrawnOnce(v, blocksOf(PeachToGrape, ItemsPitch, 7));
    await bothThemes(v);
  });

  it('one side zoomed in, the other way: grape in chr2’s 300 kb on items, peach at the fit on bins, each block once', async () => {
    const v = await open('across-bins');
    await openZoomedGrape(v);
    await vi.advanceTimersByTimeAsync(4000);
    v.draw();
    const g = v.view.get(GrapeSection)!;
    expect(g.end - g.start).toBeLessThan(400_000);
    expect(v.view.get(PeachSection)!.start).toBe(0);
    const lay = bridges(v)[0].lays[0].lay as unknown as { reads: { level?: number }[] };
    // grape the near side (items), peach's bins the rest
    expect(lay.reads.map((r) => r.level !== undefined)).toEqual([true, false]);
    expectDrawnOnce(v, blocksOf(PeachToGrape, BinsPitch, 7));
    await bothThemes(v);
  });

  it('a stub: grape without chr4 (its `include`), peach Pp03’s tail drawn as stubs, nothing unplaced', async () => {
    const v = await open('across-stub');
    const lay = bridges(v)[0].lays[0].lay as unknown as { n: number; flags: Uint8Array };
    let stubs = 0;
    for (let i = 0; i < lay.n; i++) if (lay.flags[i] & 2) stubs++;
    expect(stubs).toBeGreaterThan(5);
    expect(v.layoutIssues.value.map((i) => i.code)).not.toContain('seq.render.bridge.unplaced');
    await bothThemes(v);
  });

  it('pending, a side at a time (EC-209): the hatch on grape’s edge until its reads arrive, then ribbons', async () => {
    const v = await open('across-pending', DelayMs + 100);
    const lay = () => bridges(v)[0].lays[0].lay as unknown as { n: number; hatch: number[] };
    // peach's reads in, grape's not: its edge hatched
    expect(lay().hatch.length).toBeGreaterThan(0);
    await vi.advanceTimersByTimeAsync(30_000);
    v.draw();
    expect(lay().hatch.length).toBe(0);
    expect(lay().n).toBeGreaterThan(80);
    await bothThemes(v);
  });

  it('a reader without bridges: the band gone, the sections stacked with the sections’ gap', async () => {
    const v = await open('across-hidden');
    expect(v.frame!.layout.bridges).toBeUndefined();
    const [a, b] = v.frame!.layout.sections;
    expect(b.y - (a.y + a.height)).toBeLessThan(40);
    await bothThemes(v, 0);
  });

  it('one genome against itself (`mirror`): each duplication drawn both ways, the lower section following', async () => {
    const v = await open('across-self');
    expect(v.view.isFollowing(PeachBelowSection)).toBe(true);
    const ribbons = bridges(v)[0].ribbons();
    // Pp01's stretch on Pp03 and back
    const up = ribbons.filter((r) => named(v, r.aSystem) === 'Pp01' && named(v, r.bSystem) === 'Pp03');
    const down = ribbons.filter((r) => named(v, r.aSystem) === 'Pp03' && named(v, r.bSystem) === 'Pp01');
    expect(up.length).toBeGreaterThan(0);
    expect(down.length).toBe(up.length);
    await bothThemes(v);
  });

  it('three genomes chained (EC-235): a bridge between each two, both drawn', async () => {
    const v = await open('across-chain');
    const list = bridges(v);
    expect(list.map((b) => [b.a, b.b])).toEqual([
      [PeachSection, GrapeSection],
      [GrapeSection, PlumSection],
    ]);
    for (const b of list) expect(b.count).toBeGreaterThan(100);
    // Shift+A at grape: both neighbours brought to what grape shows, one jump
    const before = [v.view.get(PeachSection)!.end, v.view.get(PlumSection)!.end];
    v.view.zoom(GrapeSection, 8);
    await vi.advanceTimersByTimeAsync(30_000);
    v.draw();
    const jumps = v.view.jumps.value;
    expect(v.runCommand('bridge.align', { section: GrapeSection })).toBe('Aligned');
    expect(v.view.jumps.value).toBe(jumps + 1);
    expect(v.view.get(PeachSection)!.end).toBeLessThan(before[0]);
    expect(v.view.get(PlumSection)!.end).toBeLessThan(before[1]);
    await bothThemes(v, 2);
  });
});
