/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The Complex's stories, headless on a fake clock, composed as the stories compose them (a viewer of the view on the
// host's stack, its SVG export, the readout panel), over the fixtures read through the Node assets adapter: the grid
// over the complex (a hovered cell named by both residues and chains, the crosshair across the strips, a rect dragged
// across the chain boundary split at it, its summary), each editing view (one patch a step, undone and redone), the
// ligands' arcs and refs, the confidence read from AlphaFold DB's files; every view in both themes, canvas and SVG.

import { readFileSync } from 'node:fs';
import { DarkTokens, type GridLay, type PairsLay, recordDrawLists } from '@molstar/seqx-core/render';
import { ComplexViews } from '@molstar/seqx-examples';
import { copiesScale } from '@molstar/seqx-examples/complex/blocks';
import { Data, entryLigands, Files } from '@molstar/seqx-examples/complex/data';
import { nodeFixtureAssets, PackageFolder } from '@molstar/seqx-examples/node-assets';
import { Describer } from '@molstar/seqx-viewer/describe';
import type { SeqViewer } from '@molstar/seqx-viewer/viewer';
import { DomHost } from '@molstar/ui/dom';
import { afterEach, beforeEach, describe, expect, it, vi } from 'vitest';
import { composeEditing, type EditingSession, historyRows } from '../../features/editing/compose.js';
import {
  bookmarkInterface,
  hideLigandsPast,
  holdReadout,
  readout,
  setCopies,
  setPlddtBins,
  shownLigands,
} from '../complex/compose.js';

const view = (id: string) => ComplexViews.find((v) => v.id === id)!;

/** Cyclin A2's letters as its UniProt entry has them. */
const P20248 = JSON.parse(readFileSync(`${PackageFolder}${Files.cyclin}`, 'utf8')).sequence.value as string;

let session: EditingSession | undefined;
const raf = DomHost.requestAnimationFrame;

beforeEach(() => {
  DomHost.requestAnimationFrame = (fn) => {
    setTimeout(() => fn(performance.now()), 16);
  };
  vi.useFakeTimers({ toFake: ['setTimeout', 'clearTimeout', 'Date', 'performance'] });
});

afterEach(() => {
  DomHost.requestAnimationFrame = raf;
  session?.dispose();
  session = undefined;
  vi.useRealTimers();
});

/** A story's session: the view on the host's stack, its fixtures read from the package's folder. */
async function open(v = view('blocks'), build = v.build) {
  const s = composeEditing({ ...v, build }, { ports: { assets: () => nodeFixtureAssets() }, delayMs: 0 });
  session = s;
  s.viewer.resize(720, 900);
  await s.viewer.load(build());
  await vi.advanceTimersByTimeAsync(600);
  s.viewer.draw();
  return s;
}

async function svgOf(v: SeqViewer, prefix = 'cx-') {
  const p = v.renderSvg({ idPrefix: prefix });
  await vi.advanceTimersByTimeAsync(400);
  return p;
}

const centre = (a: { cell(i: number): void; sx0: number; sx1: number }, i: number) => {
  a.cell(i);
  return (a.sx0 + a.sx1) / 2;
};

describe('the grid over the complex, headless', () => {
  it('a hovered cell: both residues and chains in the hover card, the crosshair across its strips, the cell outlined once', async () => {
    const s = await open(view('complex-grid'));
    const v = s.viewer;
    // the panel holds the letters it names (read from the files): as the story mounts it
    holdReadout(v);
    await vi.advanceTimersByTimeAsync(200);
    expect(v.frame!.layout.issues).toEqual([]);
    const st = v.frame!.layout.sections[1].placement.strips;
    expect(st.n).toBe(4);
    // CDK2's K33 (x) against cyclin A2's 300 (y): strip 1
    const x = centre(st.axes[1], 32);
    const y = centre(st.ys[1]!, 299);
    v.input.move(x, y);
    v.draw();
    const h = v.hover.current;
    expect([h.system, h.position, h.ySystem, h.yPosition]).toEqual([
      '/data/cdk2-protein',
      32,
      '/data/cyclin-protein',
      299,
    ]);
    // the card (the viewer's: the panel once said it): both positions, both chains, the value
    const card = new Describer({
      spec: v.spec,
      store: v.store,
      addresses: v.addresses,
      view: v.view,
      layout: () => v.frame?.layout,
    }).describe(v.hover.current);
    expect(card?.rows).toContainEqual({ label: 'Chain x', value: 'CDK2' });
    expect(card?.rows).toContainEqual({ label: 'Chain y', value: 'cyclin A2' });
    // the value in the unit the PAE's theme declares (Å), as the panel once wrote it
    expect(card?.rows.find((r) => r.label === 'Value')?.value).toMatch(/^\d+(\.\d+)? Å$/);
    expect(card?.text).toBe(`pae, CDK2 K33 × cyclin A2 ${P20248[299]}300, value 9 Å, chain x CDK2, chain y cyclin A2`);
    expect(card?.rows).toContainEqual({ label: 'y', value: `cyclin A2 ${P20248[299]}300` });
    expect(card?.rows).toContainEqual({ label: 'x', value: 'CDK2 K33' });
    const rects = recordDrawLists(v.frame!.lists, ['hover'])
      .split('\n')
      .filter((l) => l.startsWith('rect '));
    // exactly: K33's column on the two strips of CDK2's x (one a y block), cyclin A2's row on the two strips of its y,
    // the cell outlined on its one strip, and both residues outlined on the blocks above
    const box = v.frame!.layout.sections[1].placement.box;
    const shape = (l: string) => {
      const [x, y, w, h] = l.split(' ').slice(1, 5).map(Number);
      return { x, y, w, h, fill: !l.split(' ')[5].startsWith('none'), inPlot: y >= box.y && y < box.y + box.height };
    };
    const shapes = rects.map(shape);
    const plot = shapes.filter((r) => r.inPlot);
    expect(plot.filter((r) => r.fill && r.h > 20).length).toBe(2);
    expect(plot.filter((r) => r.fill && r.w > 20).length).toBe(2);
    expect(plot.filter((r) => !r.fill).length).toBe(1);
    expect(shapes.filter((r) => !r.inPlot).length).toBe(2);
  });

  it('a rect dragged across the chain boundary is split at it: a rect a pair of chains, its summary in the panel', {
    tags: ['medium'],
  }, async () => {
    const s = await open(view('complex-grid'));
    const v = s.viewer;
    const st = v.frame!.layout.sections[1].placement.strips;
    // from CDK2 × CDK2 (280, 280) to cyclin A2 × cyclin A2 (20, 20): strips 0 and 3
    const x0 = centre(st.axes[0], 280);
    const y0 = centre(st.ys[0]!, 280);
    const x1 = centre(st.axes[3], 20);
    const y1 = centre(st.ys[3]!, 20);
    v.input.down(x0, y0);
    v.input.move((x0 + x1) / 2, (y0 + y1) / 2);
    v.input.move(x1, y1);
    v.input.up();
    const loc = v.selection.value.value.location;
    const parts = loc.rects.map((r) => `${r.a.system}:${[...r.a.ranges]} × ${r.b.system}:${[...r.b.ranges]}`);
    expect(parts.sort()).toEqual(
      [
        '/data/cdk2-protein:280,298 × /data/cdk2-protein:280,298',
        '/data/cdk2-protein:280,298 × /data/cyclin-protein:0,21',
        '/data/cyclin-protein:0,21 × /data/cdk2-protein:280,298',
        '/data/cyclin-protein:0,21 × /data/cyclin-protein:0,21',
      ].sort(),
    );
    expect(readout(v).rect).toMatch(
      /^Aligned CDK2 281–298, cyclin A2 1–21, Scored CDK2 281–298, cyclin A2 1–21, mean \d+\.\d Å$/,
    );
    // the selection drawn on each strip it covers
    v.draw();
    expect(v.frame!.lists.passes.selection.rects.n).toBeGreaterThanOrEqual(4);
  });

  it('every strip’s grid drawn in both themes, canvas and SVG alike; the boundaries hairlines', {
    tags: ['medium'],
  }, async () => {
    const s = await open(view('complex-grid'));
    const v = s.viewer;
    for (const dark of [false, true]) {
      if (dark) v.setTokens(DarkTokens);
      await vi.advanceTimersByTimeAsync(100);
      v.draw();
      const g = v.frame!.layout.sections[1];
      for (let i = 0; i < 4; i++) expect(g.lays[i].filter((l) => (l.lay as GridLay).grid).length).toBe(1);
      const r = await svgOf(v);
      expect(r.issues).toEqual([]);
      // the hairlines' sizes (the export's page may place the plot at another top)
      const hair = (lists: Parameters<typeof recordDrawLists>[0]) =>
        recordDrawLists(lists, ['labels'])
          .split('\n')
          .filter((l) => /^rect \S+ \S+ (1 \d{3,}|\d{3,} 1) /.test(l))
          .map((l) => l.split(' ').slice(3, 5).join(' '));
      expect(hair(r.lists)).toEqual(hair(v.frame!.lists));
      expect(hair(r.lists).length).toBe(2);
    }
  });
});

describe('the Complex editing views, headless', () => {
  it('blocks: the interface bookmarked, the copies split and merged: three steps, undone and redone', async () => {
    const s = await open(view('blocks'));
    const v = s.viewer;
    const labels = () => v.frame!.layout.sections[0].section.blocks!.map((b) => b.label);
    expect(labels()).toEqual(['CDK2 (A) ×2', 'cyclin A (B) ×2', 'ATP ×2']);
    expect(await bookmarkInterface(s)).not.toBe('Already bookmarked');
    expect(v.addresses.refs.value.map((r) => r.id)).toContain('interface');
    await setCopies(s, 'complex', 'each');
    await vi.advanceTimersByTimeAsync(200);
    v.draw();
    expect(labels()).toEqual(['CDK2 (A)', 'cyclin A (B)', 'CDK2 (C)', 'cyclin A (D)', 'ATP', 'ATP']);
    expect(historyRows(s).done.length).toBe(2);
    await s.stack.undo();
    await s.stack.undo();
    await vi.advanceTimersByTimeAsync(200);
    v.draw();
    expect(labels()).toEqual(['CDK2 (A) ×2', 'cyclin A (B) ×2', 'ATP ×2']);
    expect(v.addresses.refs.value.map((r) => r.id)).not.toContain('interface');
    await s.stack.redo();
    expect(v.addresses.refs.value.map((r) => r.id)).toContain('interface');
  });

  it('blocks: a hover on the merged CDK2 block lights the position in every member (one Location by system)', async () => {
    const s = await open(view('blocks'));
    const v = s.viewer;
    const sec = v.frame!.layout.sections[0];
    const letters = sec.tracks.find((t) => t.track.address.endsWith('/letters'))!.box;
    const a = sec.placement.strips.axes[0];
    v.input.move(centre(a, 44), letters.y + letters.height / 2);
    v.draw();
    expect([v.hover.current.system, v.hover.current.position]).toEqual(['/data/cdk2', 44]);
    expect(sec.section.blocks![0].members).toEqual(['/views/main/complex/a', '/views/main/complex/c']);
  });

  it('confidence: AlphaFold DB’s values drawn; other bins one theme patch, undone; cyclin A2’s first 172 gray', async () => {
    const s = await open(view('confidence'));
    const v = s.viewer;
    expect(v.frame!.layout.issues).toEqual([]);
    const lay = (section: number, track: string) =>
      v.frame!.layout.sections[section].lays[
        v.frame!.layout.sections[section].tracks.findIndex((t) => t.track.address.endsWith(`/${track}`))
      ];
    expect(lay(0, 'plddt').length).toBe(1);
    const colorOf = () =>
      JSON.stringify(
        (
          v.spec.shown.value
            .children!.find((c) => c.id === 'views')!
            .children![0].children![0].children!.find((t) => t.id === 'plddt')!.children![0].params as Record<
            string,
            unknown
          >
        ).color,
      );
    const before = colorOf();
    expect(await setPlddtBins(s, [0, 60, 80, 90, 100])).toBe('Set the confidence bins');
    expect(colorOf()).toContain('60');
    await s.stack.undo();
    expect(colorOf()).toBe(before);
    // the gray range on cyclin A2 over 1–172: an item of the unobserved layer, on its system
    const hit = v.frame!.layout.sections[1];
    const held = hit.tracks.findIndex((t) => t.track.address.endsWith('/held'));
    const box = hit.tracks[held].box;
    const a = hit.axis!;
    const owners = new Set<string>();
    const x = centre(a, 100);
    for (let y = box.y; y < box.y + box.height; y += 1) {
      const h = v.hit(x, y);
      if (h?.element) owners.add(`${h.system} ${h.element.owner}`);
    }
    expect([...owners]).toEqual(['/data/cyclin-protein /data/cyclin-unobserved']);
    const r = await svgOf(v);
    expect(r.issues).toEqual([]);
  });

  it('ligands: a chip a ligand with its arcs; an arc hovered lights both ends; the past-2 Å ones hidden by one patch, undone', async () => {
    const s = await open(view('ligands'));
    const v = s.viewer;
    const sec = v.frame!.layout.sections[0];
    expect(sec.section.blocks!.map((b) => b.label)).toEqual([
      'CDK2 (1FIN A)',
      'ATP',
      ...entryLigands().map((l) => l.comp),
    ]);
    const arcs = sec.sectionLays.map((l) => l.lay as PairsLay);
    const expected = [Data.ligands[0], ...entryLigands()].map((l) => new Set(l.contacts.map(([r]) => r)).size);
    expect(arcs.map((l) => l.n)).toEqual(expected);
    // each ligand's arcs in a color of its own (the palette's token by its place)
    const colors = new Set(
      recordDrawLists(v.frame!.lists, ['data'])
        .split('\n')
        .filter((l) => l.startsWith('path '))
        .map((l) => /#[0-9a-f]{8}/.exec(l)?.[0]),
    );
    expect(colors.size).toBe(expected.length);
    // the refs the builder wrote: one a ligand
    expect(v.addresses.refs.value.length).toBe(expected.length);
    const all = shownLigands(s);
    expect(all.length).toBe(expected.length);
    expect(await hideLigandsPast(s, 2)).toBe('Hide ligands past 2 Å');
    const past = [Data.ligands[0], ...entryLigands()]
      .filter((l) => l.resolution > 2)
      .map((l) => `lig-${l.comp.toLowerCase()}`);
    expect(past.length).toBeGreaterThan(0);
    expect(shownLigands(s)).toEqual(all.filter((id) => !past.includes(id)));
    await s.stack.undo();
    expect(shownLigands(s)).toEqual(all);
    const r = await svgOf(v);
    expect(r.issues).toEqual([]);
  });

  it('the 200 copies: one block ×200, split by one patch into 200, merged again by undo', async () => {
    const scale = { ...view('blocks'), id: 'copies-scale', build: () => copiesScale() };
    const s = await open(scale, scale.build);
    const v = s.viewer;
    const n = () => v.frame!.layout.sections[0].section.blocks!.length;
    expect(n()).toBe(1);
    expect(v.frame!.layout.sections[0].section.blocks![0].label).toBe('CDK2 1 ×200');
    await setCopies(s, 'copies', 'each');
    await vi.advanceTimersByTimeAsync(200);
    v.draw();
    expect(n()).toBe(200);
    await s.stack.undo();
    await vi.advanceTimersByTimeAsync(200);
    v.draw();
    expect(n()).toBe(1);
  });
});

// the views that take long to draw are the medium lane's; the others stay in the quick lane as its representatives
const SlowViews = new Set(['pae', 'complex-grid', 'blocks-stacked', 'blocks-stacked-split', 'blocks-ring']);

describe('every Complex view in both themes, headless', () => {
  for (const ex of ComplexViews)
    it(
      `${ex.id}: laid out and exported with no issue, light and dark`,
      SlowViews.has(ex.id) ? { tags: ['medium'] } : {},
      async () => {
        const s = await open(ex);
        for (const dark of [false, true]) {
          if (dark) s.viewer.setTokens(DarkTokens);
          await vi.advanceTimersByTimeAsync(100);
          s.viewer.draw();
          expect(s.viewer.frame!.layout.issues).toEqual([]);
          const r = await svgOf(s.viewer);
          expect(r.issues).toEqual([]);
          expect(r.settled).toBe(true);
        }
      },
    );
});
