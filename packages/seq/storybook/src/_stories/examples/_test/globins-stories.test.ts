/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The Globins stories, headless on a fake clock, composed as the stories compose them (a viewer of the view and its SVG
// export): each view loads and draws its rows (the figure its text), and its export has no issue.

import {
  type BlockPlacement,
  DarkTokens,
  type PairsLay,
  type RangesLay,
  type RowsLay,
  recordDrawLists,
  SlimRailPx,
} from '@molstar/seqx-core/render';
import { FeatureViews, GlobinsViews } from '@molstar/seqx-examples';
import { readFileSync } from 'node:fs';
import { readStockholm } from '@molstar/seqx-core/data/alignment-readers';
import { ColumnsFile, ReferenceRow, SeedFile } from '@molstar/seqx-examples/globins/data';
import { HemeHistidines, HemeUnitContacts } from '@molstar/seqx-examples/globins/hemes';
import { PackageFolder } from '@molstar/seqx-examples/node-assets';
import { LinksHoverChannel, LinksHoverProvider } from '@molstar/seqx-viewer/behaviors';
import { viewerItem } from '@molstar/seqx-viewer/config';
import type { SeqViewer } from '@molstar/seqx-viewer/viewer';
import { afterEach, describe, expect, it, vi } from 'vitest';
import { openExample } from './draw-checks.js';

let viewer: SeqViewer | undefined;

afterEach(() => {
  viewer?.retire();
  viewer = undefined;
  vi.useRealTimers();
});

function rowsLay(v: SeqViewer): RowsLay | undefined {
  for (const s of v.frame!.layout.sections)
    for (const b of s.lays) for (const l of b ?? []) if (l.type.name === 'rows') return l.lay as RowsLay;
  return undefined;
}

describe('Globins stories, headless', () => {
  it('the family: the positions rail beside the vertical overview, under the ruler, and under the last track in the bottom story', async () => {
    const family = GlobinsViews.find((v) => v.id === 'family')!;
    viewer = await openExample(family, {});
    const at = (v: SeqViewer) => v.frame!.layout.sections[0].rail!;
    const ruler = viewer.frame!.layout.sections[0].tracks[0].box;
    expect(at(viewer).bottom).toBe(ruler.y + ruler.height + SlimRailPx);
    viewer.retire();
    viewer = await openExample(family, { config: [[viewerItem('navigator.railPlace'), 'bottom']] as never });
    const last = viewer.frame!.layout.sections[0].tracks.at(-1)!.box;
    expect(at(viewer).top).toBeGreaterThanOrEqual(last.y + last.height);
    viewer.retire();
    viewer = await openExample(family, { config: [[viewerItem('navigator.railPlace'), ['ruler', 'pinned']]] as never });
    const rails = viewer.frame!.layout.sections[0].rails!;
    expect(rails.length).toBe(2);
    expect(rails[1].pinned).toBe(true);
    expect(rails[1].bottom).toBe(viewer.height);
  });

  it('hemoglobin with `links.hover` on, as its story: a heme’s chip lights the residues it contacts, its proximal histidine lights the heme', async () => {
    viewer = await openExample(GlobinsViews.find((v) => v.id === 'hemoglobin')!, {
      fixtures: true,
      behaviors: [{ provider: LinksHoverProvider, params: undefined }],
    });
    const v = viewer;
    const owner = {};
    const section = '/views/main/hb';
    const lit = () => v.highlight.get(LinksHoverChannel)?.location.parts.map((p) => [p.system, Array.from(p.ranges)]);
    v.hover.set({ section, system: '/data/hb/hem-a', position: 0 }, owner);
    // the contacts are the complex's links: each residue a run of one, 0-based
    const a = HemeUnitContacts.find((c) => c.heme === 'A')!.residues;
    const runs = lit()!;
    expect(runs.map((r) => r[0])).toEqual(['/data/hb/A']);
    const positions = new Set<number>();
    const r = runs[0][1] as number[];
    for (let k = 0; k < r.length; k += 2) for (let i = r[k]; i < r[k + 1]; i++) positions.add(i);
    expect([...positions].sort((x, y) => x - y)).toEqual(a.map((n) => n - 1));
    // the proximal histidine (His87 of α): the heme it holds
    v.hover.set({ section, system: '/data/hb/A', position: HemeHistidines.hba.proximal - 1 }, owner);
    expect(lit()).toEqual([['/data/hb/hem-a', [0, 1]]]);
    v.hover.clear(owner);
  });

  it('hemoglobin: the chains merged ×2, the hemes one chip ×4; eight helices a molecule, every heme contact an arc (per pair of units); both themes', async () => {
    viewer = await openExample(GlobinsViews.find((v) => v.id === 'hemoglobin')!, { fixtures: true });
    const v = viewer;
    expect(v.layoutIssues.value).toEqual([]);
    const sec = v.frame!.layout.sections[0];
    const p = sec.placement as BlockPlacement;
    expect(p.labels).toEqual(['α (A) ×2', 'β (B) ×2', 'heme ×4']);
    expect(p.strips.axes.map((a, i) => [a.system, p.strips.chip[i]])).toEqual([
      ['/data/hb/A', 0],
      ['/data/hb/B', 0],
      ['/data/hb/hem-a', 1],
    ]);
    // every block's lays (a block a strip)
    const lays = (track: string, type: string) =>
      sec.lays
        .flatMap((b) => b ?? [])
        .filter((l) => l.type.name === type && (l.lay as { key: string }).key.includes(`/${track}/`));
    const helices = lays('helices', 'ranges').flatMap((l) => (l.lay as RangesLay).keys);
    expect([...new Set(helices)].sort()).toEqual(['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H']);
    expect(helices.length).toBe(16);
    // a pairs layer across blocks is the section's (its arcs between the blocks' strips)
    const arcs = sec.sectionLays.reduce((n, l) => n + (l.type.name === 'pairs' ? (l.lay as PairsLay).n : 0), 0);
    expect(arcs).toBe(HemeUnitContacts.reduce((n, c) => n + c.residues.length, 0));
    // the bookmarks the builder wrote: each molecule's heme pocket
    expect(v.addresses.refs.value.map((r) => r.id)).toEqual(['pocket-a', 'pocket-b']);
    for (const tokens of [undefined, DarkTokens]) {
      if (tokens) {
        v.setTokens(tokens);
        v.draw();
      }
      const svg = v.renderSvg({ idPrefix: 'h-' });
      await vi.advanceTimersByTimeAsync(400);
      const r = await svg;
      expect(r.issues.filter((i) => i.severity === 'error')).toEqual([]);
      for (const t of ['α (A) ×2', 'β (B) ×2', 'heme ×4']) expect(r.svg).toContain(t);
    }
  });

  // `scale`'s 10,000 rows (its tiles, its export) are slow by nature: the medium lane; `lazy`'s are the quick representative
  for (const id of ['family', 'folded', 'figure', 'scale', 'lazy'])
    it(`${id}: loads, draws its rows and exports them with no issue`, {
      tags: id === 'scale' ? ['medium'] : [],
    }, async () => {
      viewer = await openExample(GlobinsViews.find((v) => v.id === id)!, { fixtures: true });
      const v = viewer;
      const lay = rowsLay(v);
      expect(lay, id).toBeDefined();
      expect(lay!.display.count).toBeGreaterThan(id === 'folded' ? 2 : 10);
      const p = v.renderSvg({ idPrefix: 'g-' });
      await vi.advanceTimersByTimeAsync(400);
      const r = await p;
      expect(r.issues.filter((i) => i.severity === 'error')).toEqual([]);
      expect(r.svg.length).toBeGreaterThan(1000);
    });

  it('family: each row’s species on the left of its numbers, its identity to the reference row on the right', async () => {
    viewer = await openExample(GlobinsViews.find((v) => v.id === 'family')!, { fixtures: true });
    const fam = rowsLay(viewer)!;
    expect(fam.cols.titleText.slice(0, fam.cols.titleN)).toEqual(['Species', 'Identity']);
    expect(fam.cols.width.every((w) => w > 0)).toBe(true);
    // the first row shown: its species from UniProt's file joined by its key, its identity a share
    const species = new Map(
      readFileSync(`${PackageFolder}${ColumnsFile}`, 'utf8')
        .trim()
        .split('\n')
        .map((l) => l.split('\t'))
        .map((c) => [c[0], c[1]] as const),
    );
    const first = fam.first;
    const cell = (k: number) => {
      for (let q = 0; q < fam.cols.n; q++)
        if (fam.cols.display[q] === first && fam.cols.col[q] === k) return fam.cols.cellText(q);
      return undefined;
    };
    const label = fam.label(fam.display.rowAt(first));
    // elided to the column's width
    expect(species.get(label)!.startsWith(cell(0)!.replace(/…$/, ''))).toBe(true);
    // the identity to the reference row as the measure says it: the same residues over the columns where either has
    // one (gap against gap and an insert left out), a percent to one decimal
    const seed = readStockholm(readFileSync(`${PackageFolder}${SeedFile}`, 'utf8')).rows;
    const ref = seed.find((r) => r.label === ReferenceRow)!.letters;
    const row = seed.find((r) => r.label === label)!.letters;
    let same = 0;
    let over = 0;
    for (let j = 0; j < row.length; j++) {
      const [a, b] = [row[j], ref[j]];
      if (a === '.' || b === '.' || a !== a.toUpperCase() || b !== b.toUpperCase()) continue;
      if (a === '-' && b === '-') continue;
      over++;
      if (a === b && a !== '-') same++;
    }
    expect(cell(1)).toBe(`${Math.round((same / over) * 1000) / 10}%`);
    // left between the label and the first numbers, right past the data
    expect(fam.cols.x[0]).toBeGreaterThan(fam.labelX);
    expect(fam.cols.x[1]).toBeGreaterThan(fam.box.x + fam.box.width);
  });

  it('family and folded: labels read whole beside the tree (its own width), and the wedges stop before them', async () => {
    viewer = await openExample(GlobinsViews.find((v) => v.id === 'family')!, { fixtures: true });
    const fam = rowsLay(viewer)!;
    expect(fam.treeX1 - fam.treeX0).toBe(120);
    expect(fam.labelText).toContain('GLB1_GLYDI/28-141');
    viewer.retire();
    viewer = await openExample(GlobinsViews.find((v) => v.id === 'folded')!, { fixtures: true });
    const folded = rowsLay(viewer)!;
    // two clades folded, each a summary row: its label "first +n", the first label cut and never its count (EC-100)
    expect(folded.foldIds.length).toBe(2);
    const summaries: string[] = [];
    for (let i = folded.first; i < folded.last; i++)
      if (folded.display.rowAt(i) < 0) summaries.push(folded.labelText[i - folded.first]);
    expect(summaries.length).toBe(2);
    for (const t of summaries) expect(t).toMatch(/\+\d+$/);
    for (let k = 0; k + 4 < folded.wedges.length; k += 5) expect(folded.wedges[k + 2]).toBeLessThan(folded.labelX);
  });

  it('scale at the fit: the rows are image tiles with their pixels drawn, and the export carries them as images', {
    tags: ['medium'],
  }, async () => {
    viewer = await openExample(GlobinsViews.find((v) => v.id === 'scale')!, { fixtures: true });
    // every row in the band: 10,000 rows in 240 px
    const key = rowsLay(viewer)!.key;
    viewer.view.setBand(key.slice(0, key.lastIndexOf('/')), { bandTop: 0, rowScale: 240 / 10_000 });
    await vi.advanceTimersByTimeAsync(50);
    viewer.draw();
    const lay = rowsLay(viewer)!;
    expect(lay.display.count * lay.rowScale).toBeLessThanOrEqual(lay.box.height + 1e-6);
    expect(lay.images.length).toBeGreaterThan(0);
    // every tile has opaque pixels (the base level holds no shares: a sample there is whole)
    for (const img of lay.images) {
      let alpha = 0;
      for (let i = 3; i < img.data.length; i += 4) alpha += img.data[i];
      expect(alpha).toBeGreaterThan(0);
    }
    const p = viewer.renderSvg({ idPrefix: 's-' });
    await vi.advanceTimersByTimeAsync(400);
    const svg = (await p).svg;
    expect(svg.match(/<image [^>]*href="data:image\/png;base64,/g)?.length).toBe(lay.images.length);
  });

  for (const id of ['rows-rail'])
    it(`${id}: the rows rail at the band's right edge, live (its rows scroll), drawn; a drag on it scrolls them`, async () => {
      viewer = await openExample(FeatureViews.find((v) => v.id === id)!, { fixtures: true });
      const v = viewer;
      const s = v.frame!.layout.sections.find((x) => x.rowsRails?.length)!;
      expect(s, id).toBeDefined();
      const r = s.rowsRails![0];
      const lay = rowsLay(v)!;
      expect(r.vertical).toBe(true);
      expect(r.active).toBe(true);
      // beside the rows' last numbers, within the canvas
      expect(r.top).toBeGreaterThanOrEqual(lay.box.x + lay.box.width);
      expect(r.bottom).toBeLessThanOrEqual(v.width);
      // drawn in the window pass: a rect within the rail's column
      const rects = v.frame!.lists.passes.window.rects;
      let inside = 0;
      for (let i = 0; i < rects.n; i++) if (rects.x[i] >= r.top && rects.x[i] + rects.w[i] <= r.bottom) inside++;
      expect(inside).toBeGreaterThan(1);
      // a drag on its window scrolls the rows
      const c = v.context.renderConfig.navigator.window;
      const y = r.x + Math.max(c.minPx, ((r.end - r.start) / r.total) * r.width) / 2;
      const x = (r.top + r.bottom) / 2;
      expect(v.input.down(x, y)).toBe(true);
      v.input.move(x, y + r.width / 4);
      v.input.up();
      expect(v.view.band(lay.key.slice(0, lay.key.lastIndexOf('/')))?.bandTop).toBeGreaterThan(r.total / 8);
    });

  for (const id of ['family', 'folded'])
    it(`${id}: the vertical overview declared, the rows rail beside it; its window has the rows and the columns`, async () => {
      viewer = await openExample(GlobinsViews.find((v) => v.id === id)!, { fixtures: true });
      const v = viewer;
      expect(v.frame!.layout.sections.some((x) => x.rowsRails?.length)).toBe(true);
      const nav = v.navigatorModel('overview')!;
      nav.draw();
      await vi.advanceTimersByTimeAsync(600);
      v.draw();
      const w = nav.draw()!.windows[0];
      expect(w.shown).toBe(true);
      // four grips: the rows (top, bottom) and the columns (left, right)
      expect(w.axes).toBe(3);
      expect(nav.width).toBe(96);
    });

  it('folded at 1280 px: the overview’s picture draws no row columns, so no "2 hidden" chip over its top', async () => {
    viewer = await openExample(GlobinsViews.find((v) => v.id === 'folded')!, { fixtures: true });
    const v = viewer;
    v.resize(1280 - 96);
    await vi.advanceTimersByTimeAsync(600);
    v.draw();
    // the canvas's rows keep their columns, every one shown at this width
    const lay = rowsLay(v)!;
    expect(lay.cols.width.length).toBeGreaterThan(0);
    expect(lay.cols.chip).toBe('');
    const nav = v.navigatorModel('overview')!;
    nav.draw();
    await vi.advanceTimersByTimeAsync(600);
    v.draw();
    // the picture (96 px): the rows summarized, no header, no columns and no chip naming the columns it would hide
    const picture = recordDrawLists(nav.draw()!.lists);
    expect(picture).not.toContain('hidden');
  });
});
