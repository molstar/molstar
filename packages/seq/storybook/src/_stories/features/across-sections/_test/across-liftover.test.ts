/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { DarkTokens, LightTokens, type RangesLay } from '@molstar/seqx-core/render';
import { loadBridges } from '@molstar/seqx-core/render/bridge-load';
import { loadGenome } from '@molstar/seqx-core/resolve';
import {
  chainRuns,
  Deleted,
  Duplicated,
  liftImages,
  NewSection,
  OldSection,
  roseGenes,
} from '@molstar/seqx-examples/across/liftover';
import { AcrossViews } from '@molstar/seqx-examples/across/views';
import { describerFor } from '@molstar/seqx-viewer/describe';
import { goToLocus } from '@molstar/seqx-viewer/locus';
import type { SeqViewer } from '@molstar/seqx-viewer/viewer';
import { afterEach, beforeAll, describe, expect, it, vi } from 'vitest';
import { openExample } from '../../../examples/_test/draw-checks.js';

// Features/Across sections/Liftover, headless on a fake clock as the story composes it (a viewer of the view, its SVG
// export light and dark): Rose v1's genes drawn on v2 through the chain read by region (EC-212), each gene's images
// its drawn parts (a duplicated gene twice, EC-20; a deleted one not at all and a half-deleted one in part, EC-21; an
// inverted one in its order), the chain's runs as ribbons with no ribbon over the deletion, and a hover on the
// duplicated gene lighting both images, its card listing them.

let viewer: SeqViewer | undefined;

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

const view = AcrossViews.find((v) => v.id === 'across-liftover')!;

async function open(): Promise<SeqViewer> {
  viewer = await openExample(view, { delayMs: 300, wait: 30_000 });
  viewer.draw();
  return viewer;
}

const named = (v: SeqViewer, system: string) => v.view.model.value.systems.get(system)?.name ?? system;

/** The genes a section drew, by chromosome: `id start end`, each item's adjacent rows joined, sorted. */
function drawn(v: SeqViewer, section: string): Record<string, string[]> {
  const sec = v.frame!.layout.sections.find((s) => s.section.address === section)!;
  const strips = sec.placement.strips;
  const bands = sec.lays.length / strips.n;
  const out: Record<string, string[]> = {};
  for (let i = 0; i < strips.n; i++) {
    const rows: [string, number, number][] = [];
    for (let k = 0; k < bands; k++)
      for (const l of sec.lays[i * bands + k] ?? []) {
        if (!l.key.endsWith('/genes/genes')) continue;
        const r = (l.lay as RangesLay).rows;
        for (let j = 0; j < r.n; j++) rows.push([String(r.keys[j]), r.starts[j], r.ends[j]]);
      }
    rows.sort((x, y) => x[0].localeCompare(y[0]) || x[1] - y[1]);
    const joined: [string, number, number][] = [];
    for (const r of rows) {
      const last = joined[joined.length - 1];
      if (last && last[0] === r[0] && last[2] >= r[1]) last[2] = Math.max(last[2], r[2]);
      else joined.push([...r]);
    }
    out[named(v, strips.axes[i].system)] = joined.map((r) => r.join(' ')).sort();
  }
  return out;
}

/** Every gene's images through the chain, by v2 chromosome: the truth, computed from the runs. */
function truth(): Record<string, string[]> {
  const runs = chainRuns();
  const out: Record<string, string[]> = { chr1: [], chr2: [], chr3: [] };
  for (const g of roseGenes())
    for (const m of liftImages(runs, g.chr, g.start, g.end)) out[m.chr].push(`${g.id} ${m.start} ${m.end}`);
  for (const k of Object.keys(out)) out[k].sort();
  return out;
}

async function bothThemes(v: SeqViewer) {
  const out: string[] = [];
  for (const dark of [false, true]) {
    const r = v.renderSvg({ idPrefix: dark ? 'ld-' : 'll-', tokens: dark ? DarkTokens : LightTokens });
    await vi.advanceTimersByTimeAsync(4000);
    const svg = await r;
    expect(svg.settled).toBe(true);
    expect(svg.issues.map((i) => `${i.code}: ${i.message}`)).toEqual([]);
    expect(svg.svg.match(/data-bridge=/g)?.length ?? 0).toBe(1);
    out.push(svg.svg);
  }
  expect(out[0]).not.toBe(out[1]);
}

describe('Features/Across sections/Liftover, headless', () => {
  it('v1’s genes on v2 are their images through the chain: twice, in part, not at all (EC-20, EC-21, EC-212)', async () => {
    const v = await open();
    const got = drawn(v, NewSection);
    expect(got).toEqual(truth());
    const all = Object.entries(got).flatMap(([chr, rows]) => rows.map((r) => `${chr} ${r}`));
    // the duplicated gene in both places, the deleted one nowhere
    expect(new Set(all.filter((r) => r.split(' ')[1] === 'dup').map((r) => r.split(' ')[0]))).toEqual(
      new Set(['chr1', 'chr3']),
    );
    expect(all.some((r) => r.split(' ')[1] === 'lost')).toBe(false);
    // the chain's runs as ribbons: none over the deletion, the duplicated stretch's to two chromosomes
    const ribbons = v.frame!.layout.bridges!.list[0].ribbons();
    const on = (r: (typeof ribbons)[number]) => [named(v, r.aSystem), named(v, r.bSystem)];
    expect(ribbons.some((r) => on(r)[0] === 'chr2' && r.aEnd > Deleted[0] && r.aStart < Deleted[1])).toBe(false);
    const dup = ribbons.filter((r) => on(r)[0] === 'chr1' && r.aStart >= Duplicated[0] && r.aEnd <= Duplicated[1]);
    expect(new Set(dup.map((r) => on(r)[1]))).toEqual(new Set(['chr1', 'chr3']));
    await bothThemes(v);
  });

  it('a hover on the duplicated gene lights both images; its card lists them ("2 places")', async () => {
    const v = await open();
    await goToLocus(v, 'chr1:1,480,001-1,560,000', { section: OldSection });
    await vi.advanceTimersByTimeAsync(4000);
    v.draw();
    const sec = v.frame!.layout.sections.find((s) => s.section.address === OldSection)!;
    let hit: ReturnType<SeqViewer['hit']>;
    for (let y = sec.y; y < sec.y + sec.height && !hit; y += 2)
      for (let x = 0; x < 720 && !hit; x += 4) {
        const h = v.hit(x, y);
        if (h?.element?.key === 'dup') hit = h;
      }
    expect(hit).toBeDefined();
    v.hover.set({ ...hit!, section: undefined } as never, 'test');
    const lit = (v.hover.projectedAll()?.parts ?? []).filter((p) => p.system.startsWith('/data/v2/'));
    expect(new Set(lit.map((p) => named(v, p.system)))).toEqual(new Set(['chr1', 'chr3']));
    const card = describerFor({ ...v, layout: () => v.frame?.layout }).describe(v.hover.current)!;
    const places = card.rows.find((r) => r.label === 'Places');
    expect(places?.value).toMatch(/^2 places: chr1 .*, chr3 /);
  });

  it('the inverted gene reads `−` on v2 and `+` on v1; its thick part lands inside its image (EC-21)', async () => {
    const v = await open();
    const strandOn = async (section: string, locus: string) => {
      await goToLocus(v, locus, { section });
      await vi.advanceTimersByTimeAsync(20_000);
      v.draw();
      v.draw();
      const sec = v.frame!.layout.sections.find((s) => s.section.address === section)!;
      let hit: ReturnType<SeqViewer['hit']>;
      for (let y = sec.y; y < sec.y + sec.height && !hit; y += 2)
        for (let x = 0; x < v.width && !hit; x += 2) {
          const h = v.hit(x, y);
          if (h?.element?.key === 'inv') hit = h;
        }
      expect(hit).toBeDefined();
      v.hover.set({ ...hit!, section: undefined } as never, 'test');
      const card = describerFor({ ...v, layout: () => v.frame?.layout }).describe(v.hover.current)!;
      return card.rows.find((r) => r.label === 'Strand')?.value;
    };
    const inv = roseGenes().find((g) => g.id === 'inv')!;
    const [image] = liftImages(chainRuns(), inv.chr, inv.start, inv.end);
    const [thick] = liftImages(chainRuns(), inv.chr, inv.thickStart, inv.thickEnd);
    expect(thick.start).toBeGreaterThan(image.start);
    expect(thick.end).toBeLessThan(image.end);
    const fmt = (n: number) => n.toLocaleString('en-US');
    expect(await strandOn(NewSection, `chr2:${fmt(image.start - 9_999)}-${fmt(image.end + 10_000)}`)).toBe('−');
    expect(await strandOn(OldSection, `chr2:${fmt(inv.start - 9_999)}-${fmt(inv.end + 10_000)}`)).toBe('+');
    // the lifted row's thick part as drawn: the image's own, through the inverted run
    const sec = v.frame!.layout.sections.find((s) => s.section.address === NewSection)!;
    const lifted = sec.lays
      .flat()
      .find((l) => l?.key.endsWith('/genes/genes') && (l.lay as RangesLay).rows.keys.includes('inv'));
    const rows = (lifted!.lay as RangesLay).rows;
    const j = rows.keys.indexOf('inv');
    const col = (name: string) => rows.fields[j]?.find((c) => c.name === name)?.column as ArrayLike<number>;
    expect([col('thick_start')[rows.row[j]], col('thick_end')[rows.row[j]]]).toEqual([thick.start, thick.end]);
  });
});
