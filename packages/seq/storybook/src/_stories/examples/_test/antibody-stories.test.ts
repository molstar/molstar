/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { readFileSync } from 'node:fs';
import { hashNode, resolveAddress } from '@molstar/core/repr/tree';
import { DarkTokens, LightTokens, type RowsLay, type RulerLay, recordDrawLists } from '@molstar/seqx-core/render';
import { AntibodyViews } from '@molstar/seqx-examples';
import { theDesign, wildType } from '@molstar/seqx-examples/antibody/data';
import { PackageFolder } from '@molstar/seqx-examples/node-assets';
import { DomHost } from '@molstar/ui/dom';
import { afterEach, beforeEach, describe, expect, it, vi } from 'vitest';
import { composeEditing, type EditingSession } from '../../features/editing/compose.js';
import {
  bookmark,
  bookmarkLines,
  bookmarks,
  labelsAt,
  readoutLines,
  recolor,
  remove,
  rename,
  selectCdr,
  setRuler,
} from '../antibody/compose.js';

// The Antibody's stories, headless on a fake clock, composed as the stories compose them (examples.md, Antibody:
// Views): the rulers' labels in both schemes drawn against the data (the wild type's Kabat insertions, the design's
// IMGT 111.1 and 112.1) in the draw lists and the SVG, both themes; the hover panel's labels in both numberings; a query
// by label, and a range written backwards in IMGT an issue (EC-7); the ruler's scheme switched, and the CDR bookmarks
// added, renamed, recolored and removed, each one step undone and redone.

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

/** The fixtures read at once (the file system would not settle on the fake clock). */
const fixtures = () => ({
  reader: {
    readText: async (uri: string) => readFileSync(`${PackageFolder}${uri}`, 'utf8'),
    readBytes: async (uri: string) => new Uint8Array(readFileSync(`${PackageFolder}${uri}`)),
  },
});

const antibody = (id: string) => AntibodyViews.find((v) => v.id === id)!;

async function open(id: string, width = 1200) {
  const view = antibody(id);
  const s = composeEditing(view, { ports: { assets: fixtures }, delayMs: 0 });
  session = s;
  // tall enough for both chains' sections
  s.viewer.resize(width, 900);
  await s.viewer.load(view.build());
  await vi.advanceTimersByTimeAsync(300);
  s.viewer.draw();
  return s;
}

const settle = async (s: EditingSession) => {
  await vi.advanceTimersByTimeAsync(100);
  s.viewer.draw();
};

const hash = (s: EditingSession) => hashNode(s.viewer.spec.root.value);

/**
 * What a story draws over a window of both chains: the canvas's draw lists' texts and the SVG's (the viewer's export,
 * over its store, as the story's SVG pane is drawn).
 */
async function drawn(id: string, window: [number, number], tokens = LightTokens) {
  const s = await open(id, 2400);
  for (const c of ['light', 'heavy']) s.viewer.view.setRange(`/views/main/${c}`, ...window);
  await settle(s);
  const f = s.viewer.draw() ?? s.viewer.frame!;
  const texts = recordDrawLists(f.lists, ['data', 'labels'])
    .split('\n')
    .filter((l) => l.startsWith('text'));
  const r = await s.viewer.renderSvg({ tokens });
  const svg = [...r.svg.matchAll(/<text[^>]*>([^<]*)<\/text>/g)].map((m) => m[1]);
  const issues = [...f.layout.issues, ...r.issues];
  session?.dispose();
  session = undefined;
  // the heavy chain's IMGT ruler: its ticks (a label, its x, whether it clears its neighbours and is drawn)
  const heavy = f.layout.sections[1].tracks.find((t) => t.track.address.endsWith('/imgt'))!;
  const lay = heavy.layers[0].lay as RulerLay;
  const ticks = Array.from({ length: lay.n }, (_, j) => ({
    label: lay.labels[j],
    x: lay.x[j],
    shown: lay.shown[j] === 1,
  }));
  return { issues, texts, svg, raw: r.svg, ticks };
}

describe('Antibody: the rulers against the data', () => {
  for (const t of ['light', 'dark'] as const)
    it(`${t}: the wild type’s Kabat ruler draws the heavy chain’s insertions, read from the file; no issue`, async () => {
      const tokens = t === 'dark' ? DarkTokens : LightTokens;
      // CDR-H2 to FR3's 82A–C, then CDR-H3: wide enough for a label a residue
      for (const [window, labels] of [
        [[46, 60], ['52A']],
        [
          [78, 92],
          ['82A', '82B', '82C'],
        ],
        [
          [98, 112],
          ['100A', '100B', '100C'],
        ],
      ] as const) {
        const r = await drawn('schemes', [...window], tokens);
        expect(r.issues).toEqual([]);
        for (const l of labels) {
          expect(r.svg).toContain(l);
          expect(r.texts.some((x) => x.includes(`"${l}"`))).toBe(true);
        }
      }
      const one = await drawn('schemes', [98, 112], tokens);
      const other = await drawn('schemes', [98, 112], t === 'dark' ? LightTokens : DarkTokens);
      expect(other.raw).not.toBe(one.raw);
    });

  it('the design’s IMGT ruler ticks 111.1 and 112.1 before 112 (inline labels, EC-7); Kabat’s 100D and 100E drawn', async () => {
    const r = await drawn('schemes-design', [100, 116]);
    expect(r.issues).toEqual([]);
    for (const l of ['111.1', '100D', '100E']) expect(r.svg).toContain(l);
    const x = (l: string) => r.ticks.find((t) => t.label === l)!.x;
    expect(x('111.1')).toBeLessThan(x('112.1'));
    expect(x('112.1')).toBeLessThan(x('112'));
    // at the closest zoom (32 px a residue) three labels a residue apart do not clear: 112.1 yields to its neighbours
    expect(r.ticks.find((t) => t.label === '112.1')!.shown).toBe(false);
  });

  it('the design against the wild type: two rows a chain, as wide as their alignment', async () => {
    const s = await open('against');
    const d = theDesign();
    const sections = s.viewer.frame!.layout.sections;
    expect(sections.map((x) => x.id)).toEqual(['light', 'heavy']);
    expect(s.viewer.frame!.layout.issues).toEqual([]);
    for (const [k, c] of (['light', 'heavy'] as const).entries()) {
      const rows = sections[k].lays.flatMap((b) => b ?? []).find((l) => l.type.name === 'rows')!.lay as RowsLay;
      expect(rows.total).toBe(2);
      expect(rows.columns).toBe(d[c].aligned.design.length);
    }
    // the gap columns where the edit list puts them: after the heavy chain's F (Kabat 100) two inserted residues, a
    // gap in the wild type's row; at the light chain's N30, deleted, a gap in the design's (no indel before either)
    const rowsOf = (c: string) =>
      (resolveAddress(s.viewer.spec.root.value, ['data', `${c}-pair`])?.params?.rows as { letters: string[] }).letters;
    const gaps = (row: string) => [...row].flatMap((x, i) => (x === '-' ? [i] : []));
    const wt = wildType();
    const f = wt.heavy.kabat.indexOf('100');
    expect(gaps(rowsOf('heavy')[1])).toEqual([f + 1, f + 2]);
    expect(gaps(rowsOf('heavy')[0])).toEqual([]);
    expect(gaps(rowsOf('light')[0])).toEqual([wt.light.kabat.indexOf('30')]);
    expect(gaps(rowsOf('light')[1])).toEqual([]);
  });
});

describe('Antibody: labels in both numberings', () => {
  it('a hovered residue: its Kabat and IMGT labels in the panel, once the files’ labels are in', async () => {
    const s = await open('schemes');
    const wt = wildType();
    const i = wt.heavy.kabat.indexOf('100A');
    expect(labelsAt(s.viewer, '/data/heavy', i)).toEqual({ kabat: '100A', imgt: '113' });
    const section = s.viewer.frame!.layout.sections[1];
    const letters = section.tracks[2];
    const x = (section.axis!.x(i) + section.axis!.x(i + 1)) / 2;
    s.viewer.input.move(x, letters.box.y + letters.box.height / 2);
    await settle(s);
    expect(readoutLines(s.viewer)[0]).toBe(`Hover: heavy ${i + 1}, Kabat 100A, IMGT 113`);
  });

  it('a query by label in either scheme selects; a range written backwards in IMGT is an issue with its hint (EC-7)', async () => {
    const s = await open('schemes-design');
    const d = theDesign();
    const sel = (a: string) => s.viewer.addresses.resolve(a);
    const k = sel('/data/heavy/(at :label 100A :numbering kabat)');
    expect(k.status).toBe('resolved');
    const i = d.heavy.kabat.indexOf('100A');
    expect(k.result?.type === 'location' && Array.from(k.result.location.parts[0].ranges)).toEqual([i, i + 1]);
    const j = d.heavy.imgt.indexOf('112.1');
    const m = sel('/data/heavy/(at :label "112.1" :numbering imgt)');
    expect(m.result?.type === 'location' && Array.from(m.result.location.parts[0].ranges)).toEqual([j, j + 1]);
    // 112 comes after 112.1 in the chain: from 112 to 112.1 runs backwards
    const back = sel('/data/heavy/(range :from 112 :to "112.1" :numbering imgt)');
    const issue = back.issues.find((x) => x.code === 'seq.query.range.backwards');
    expect(issue?.message).toBe('A range ends before it starts');
    expect(issue?.hint).toBe('write :from at or before :to, in the numbering’s order of positions');
    // in index order it is a range of two positions
    const fwd = sel('/data/heavy/(range :from "112.1" :to 112 :numbering imgt)');
    expect(fwd.result?.type === 'location' && Array.from(fwd.result.location.parts[0].ranges)).toEqual([j, j + 2]);
  });
});

describe('Antibody: the editing views', () => {
  it('the rulers’ scheme switched: one patch, one step, undone to the start', async () => {
    const s = await open('schemes');
    const start = hash(s);
    const layer = (c: string) => resolveAddress(s.viewer.spec.root.value, ['views', 'main', c, 'ruler', 'ruler']);
    expect(await setRuler(s, 'imgt', ['heavy'])).toBe('Ruler in imgt');
    expect(layer('heavy')?.params?.numbering).toBe('imgt');
    expect(layer('light')?.params?.numbering).toBe('kabat');
    // its track's label says so
    expect(resolveAddress(s.viewer.spec.root.value, ['views', 'main', 'heavy', 'ruler'])?.params?.label).toBe(
      'Heavy (B), IMGT',
    );
    await s.stack.undo();
    await settle(s);
    expect(hash(s)).toBe(start);
    await s.stack.redo();
    expect(layer('heavy')?.params?.numbering).toBe('imgt');
    // both chains in one step
    await setRuler(s, 'ordinal');
    expect([layer('light')?.params?.numbering, layer('heavy')?.params?.numbering]).toEqual(['ordinal', 'ordinal']);
    await s.stack.undo();
    expect(layer('light')?.params?.numbering).toBe('kabat');
  });

  it('a CDR bookmarked from the selection, renamed, recolored, removed; each a step, undone and redone', async () => {
    const s = await open('cdrs');
    const start = hash(s);
    expect(selectCdr(s, 'heavy', 3)).toBe('Selected');
    expect(s.viewer.addresses.selectionAddress()).toMatch(/^\/data\/heavy\/.*kabat/);
    expect(await bookmark(s, 'cdr-h3', 'CDR-H3')).toBe('Bookmarked');
    const added = hash(s);
    expect(bookmarks(s)).toEqual([{ id: 'cdr-h3', label: 'CDR-H3', color: undefined }]);
    // the bookmark names the CDR's residues: Kabat 95 to 102 of the heavy chain
    const ref = s.viewer.addresses.resolve('/refs/cdr-h3');
    const wt = wildType();
    const at = wt.heavy.kabat.indexOf('95');
    expect(ref.result?.type === 'location' && Array.from(ref.result.location.parts[0].ranges)).toEqual([
      at,
      wt.heavy.kabat.indexOf('102') + 1,
    ]);
    const refNode = () => resolveAddress(s.viewer.spec.root.value, ['refs', 'cdr-h3']);
    const to = refNode()?.params?.to;
    expect(await rename(s, 'cdr-h3', 'Heavy CDR3')).toBe('Labeled');
    // the label alone changes: what it names stays
    expect(refNode()?.params).toEqual({ ...refNode()?.params, label: 'Heavy CDR3', to });
    expect(bookmarks(s)).toEqual([{ id: 'cdr-h3', label: 'Heavy CDR3', color: undefined }]);
    await recolor(s, 'cdr-h3', '#d62728');
    expect(bookmarks(s)).toEqual([{ id: 'cdr-h3', label: 'Heavy CDR3', color: '#d62728' }]);
    // the panel lists it
    expect(bookmarkLines(s.viewer).slice(-2)).toEqual([
      'Bookmarks (1):',
      `  Heavy CDR3 [cdr-h3] #d62728: ${Array.isArray(to) ? to.join(', ') : to}`,
    ]);
    expect(await remove(s, 'cdr-h3')).toBe('Removed');
    // the ref node is gone from the spec, the selection kept
    expect(refNode()).toBeUndefined();
    expect(bookmarks(s)).toEqual([]);
    expect(bookmarkLines(s.viewer).slice(-1)).toEqual(['  (none)']);
    expect(s.viewer.addresses.selectionAddress()).toMatch(/kabat/);
    // undone one at a time: removed, recolored, renamed, added
    await s.stack.undo();
    expect(bookmarks(s)[0].color).toBe('#d62728');
    await s.stack.undo();
    expect(bookmarks(s)[0].color).toBeUndefined();
    await s.stack.undo();
    expect(bookmarks(s)[0].label).toBe('CDR-H3');
    expect(hash(s)).toBe(added);
    await s.stack.undo();
    expect(bookmarks(s)).toEqual([]);
    expect(hash(s)).toBe(start);
    await s.stack.redo();
    expect(bookmarks(s).map((b) => b.id)).toEqual(['cdr-h3']);
  });
});
