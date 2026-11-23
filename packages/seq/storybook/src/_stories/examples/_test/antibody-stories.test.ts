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
import { Describer } from '@molstar/seqx-viewer/describe';
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

// The Antibody's stories, headless on a fake clock, composed as the stories compose them: the rulers' labels in both
// schemes drawn against the data (the wild type's Kabat insertions, the design's IMGT 111.1 and 112.1) in the draw
// lists and the SVG, both themes; the hover panel's labels in both numberings; a query by label, and a range written
// backwards in IMGT an issue (EC-7); the ruler's scheme switched, and the CDR bookmarks added, renamed, recolored and
// removed, each one step undone and redone.

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
  // the ruler ticks the step's round values, and an insertion code's label where it clears them: at the closest zoom
  // (32 px a residue) CDR-H2's 52A, FR3's 82A–C and CDR-H3's 100A–C show beside the round values; at the fit only
  // where one clears them
  for (const t of ['light', 'dark'] as const)
    it(
      `${t}: the wild type’s Kabat ruler labels the heavy chain’s insertions at the closest zoom, read from the file; no issue`,
      t === 'dark' ? { tags: ['medium'] } : {},
      async () => {
        const tokens = t === 'dark' ? DarkTokens : LightTokens;
        const has = (r: { texts: string[] }, l: string) => r.texts.some((x) => x.includes(`"${l}"`));
        for (const [window, round, insertions] of [
          [[46, 60], ['50', '52'], ['52A']],
          [
            [78, 92],
            ['80', '82'],
            ['82A', '82B', '82C'],
          ],
          [
            [98, 112],
            ['98', '100'],
            ['100A', '100B', '100C'],
          ],
        ] as const) {
          const r = await drawn('schemes', [...window], tokens);
          expect(r.issues).toEqual([]);
          for (const l of [...round, ...insertions]) expect(has(r, l), l).toBe(true);
          for (const l of insertions) expect(r.svg).toContain(l);
        }
        // at the fit the round values, and an insertion only where it clears them: none crowds a round value
        const fit = await drawn('schemes', [0, 220], tokens);
        expect(has(fit, '100')).toBe(true);
        const spans = fit.ticks
          .filter((t) => t.shown)
          .map((t) => t.x)
          .sort((a, b) => a - b);
        for (let k = 1; k < spans.length; k++) expect(spans[k] - spans[k - 1]).toBeGreaterThan(4);
        const one = await drawn('schemes', [98, 112], tokens);
        const other = await drawn('schemes', [98, 112], t === 'dark' ? LightTokens : DarkTokens);
        expect(other.raw).not.toBe(one.raw);
      },
    );

  it('past the variable domain a ruler is marked, not blank: no tick, one note, a dashed line; canvas and SVG alike', async () => {
    const s = await open('schemes', 2400);
    for (const c of ['light', 'heavy']) s.viewer.view.setRange(`/views/main/${c}`, 0, 240);
    await settle(s);
    const f = s.viewer.draw() ?? s.viewer.frame!;
    const r = await s.viewer.renderSvg({ tokens: LightTokens });
    const wt = wildType();
    for (const [k, c] of (['light', 'heavy'] as const).entries()) {
      const covered = wt[c].kabat.findIndex((l) => l === null);
      expect(covered).toBeGreaterThan(90);
      for (const track of f.layout.sections[k].tracks) {
        const lay = track.layers[0].lay as RulerLay;
        if (track.layers[0].type.name !== 'ruler') continue;
        expect(lay.nruns, track.track.address).toBe(1);
        for (let j = 0; j < lay.n; j++) expect(lay.axis.positionAt(lay.x[j])).toBeLessThan(wt[c].kabat.length);
        // no tick past the last labelled position
        const last = wt[c].kabat.length - 1 - [...wt[c].kabat].reverse().findIndex((l) => l !== null);
        for (let j = 0; j < lay.n; j++) expect(lay.axis.positionAt(lay.x[j])).toBeLessThanOrEqual(last);
        expect(lay.note).toBe('not numbered');
      }
    }
    const texts = [...r.svg.matchAll(/<text[^>]*>([^<]*)<\/text>/g)].map((m) => m[1]);
    // a note for each of the four rulers, a muted dashed line for each
    expect(texts.filter((t) => t === 'not numbered')).toHaveLength(4);
    expect(r.svg.match(/stroke-dasharray="3"/g)?.length).toBeGreaterThanOrEqual(4);
    expect(f.layout.issues).toEqual([]);
  });

  it('the design’s IMGT ruler ticks 111.1 and 112.1 before 112 (inline labels, EC-7); Kabat’s 100D and 100E drawn', async () => {
    const r = await drawn('schemes-design', [100, 116]);
    expect(r.issues).toEqual([]);
    for (const l of ['112.1', '100D', '100E']) expect(r.svg).toContain(l);
    const x = (l: string) => r.ticks.find((t) => t.label === l)!.x;
    expect(x('111.1')).toBeLessThan(x('112.1'));
    expect(x('112.1')).toBeLessThan(x('112'));
    // at the closest zoom (32 px a residue) labels a residue apart do not all clear: the insertion nearest a drawn round
    // value wins (112.1 beside 112), its neighbour yields
    expect(r.ticks.find((t) => t.label === '112.1')!.shown).toBe(true);
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
  it('a hovered residue: its Kabat and IMGT labels in the hover card, once the files’ labels are in', async () => {
    const s = await open('schemes');
    const wt = wildType();
    const i = wt.heavy.kabat.indexOf('100A');
    expect(labelsAt(s.viewer, '/data/heavy', i)).toEqual({ kabat: '100A', imgt: '113' });
    const section = s.viewer.frame!.layout.sections[1];
    const letters = section.tracks[2];
    const x = (section.axis!.x(i) + section.axis!.x(i + 1)) / 2;
    s.viewer.input.move(x, letters.box.y + letters.box.height / 2);
    await settle(s);
    // the card's rows (the viewer's, not the story's): the position in each numbering of its system
    const card = new Describer({
      spec: s.viewer.spec,
      store: s.viewer.store,
      addresses: s.viewer.addresses,
      view: s.viewer.view,
      layout: () => s.viewer.frame?.layout,
    }).describe(s.viewer.hover.current);
    expect(card?.rows).toContainEqual({ label: 'kabat', value: '100A' });
    expect(card?.rows).toContainEqual({ label: 'imgt', value: '113' });
    expect(card?.title).toBe(`1N8Z_B ${i + 1}`);
    // the panel says what the card does not: the selection
    expect(readoutLines(s.viewer)).toEqual(['Selection: (none)']);
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
    // the second ruler takes the other scheme: the two never read the same
    const second = resolveAddress(s.viewer.spec.root.value, ['views', 'main', 'heavy', 'imgt']);
    expect(second?.params?.label).toBe('Heavy (B), Kabat');
    expect(
      resolveAddress(s.viewer.spec.root.value, ['views', 'main', 'heavy', 'imgt', 'ruler'])?.params?.numbering,
    ).toBe('kabat');
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

  it('the design’s rulers take the same switch (Design numberings): its ticks stay valid in IMGT, undone to Kabat', async () => {
    const s = await open('schemes-design');
    const layer = (c: string) => resolveAddress(s.viewer.spec.root.value, ['views', 'main', c, 'ruler', 'ruler']);
    expect(await setRuler(s, 'imgt')).toBe('Ruler in imgt');
    await settle(s);
    expect([layer('light')?.params?.numbering, layer('heavy')?.params?.numbering]).toEqual(['imgt', 'imgt']);
    expect(s.viewer.frame!.layout.issues).toEqual([]);
    await s.stack.undo();
    expect(layer('heavy')?.params?.numbering).toBe('kabat');
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

describe('Antibody: the rulers’ tick labels never overlap', () => {
  it('Kabat’s and IMGT’s rulers of both chains, at every zoom over the insertions: shown labels clear each other', async () => {
    const s = await open('schemes', 1200);
    const overlaps: string[] = [];
    let pairs = 0;
    for (const chain of ['light', 'heavy'])
      for (const [start, end] of [
        [0, 130],
        [60, 130],
        [70, 100],
        [76, 94],
        [78, 92],
        [80, 90],
        [98, 112],
        [100, 110],
      ] as const) {
        s.viewer.view.setRange(`/views/main/${chain}`, start, end);
        await settle(s);
        const section = s.viewer.frame!.layout.sections.find((x) => x.id === chain)!;
        for (const t of section.tracks) {
          const lay = t.layers[0]?.lay as RulerLay | undefined;
          if (!lay?.labels || t.track.address.endsWith('/letters')) continue;
          const shown = Array.from({ length: lay.n }, (_, j) => j).filter((j) => lay.shown[j] === 1);
          shown.sort((a, b) => lay.x[a] - lay.x[b]);
          for (let k = 1; k < shown.length; k++) {
            const a = shown[k - 1];
            const b = shown[k];
            const gap = lay.x[b] - lay.widths[b] / 2 - (lay.x[a] + lay.widths[a] / 2);
            pairs++;
            if (gap < 4)
              overlaps.push(
                `${chain} ${t.track.address} [${start},${end}] ${lay.labels[a]} ${lay.labels[b]} ${gap.toFixed(1)}`,
              );
          }
        }
      }
    expect(pairs).toBeGreaterThan(50);
    expect(overlaps).toEqual([]);
  });
});

describe('Antibody: the CDR labels at the fit', () => {
  // each CDR a short label (the track says the scheme): the three of a chain drawn at every fit width, none dropped
  // for crowding its neighbour (CDR-H1 by Kabat, CDR-L1 and CDR-L2 by IMGT went at narrower frames)
  it('every CDR of both chains is labelled in both schemes at a narrow and a wide frame', async () => {
    for (const width of [700, 800, 1000, 1200]) {
      const s = await open('schemes', width);
      const f = s.viewer.draw() ?? s.viewer.frame!;
      const texts = recordDrawLists(f.lists, ['data', 'labels'])
        .split('\n')
        .filter((l) => l.startsWith('text'));
      for (const chain of ['L', 'H'])
        for (const n of [1, 2, 3]) {
          const drawn = texts.filter((t) => t.startsWith(`text "CDR-${chain}${n}`)).length;
          expect(drawn, `${width}: CDR-${chain}${n} (one a scheme)`).toBe(2);
        }
      session?.dispose();
      session = undefined;
    }
  });
});

describe('Antibody: the short CDR labels keep their scheme in the card', () => {
  it('hovering CDR-H1 in the Kabat track and in the IMGT track: the cards name the track, so they differ', async () => {
    const s = await open('schemes');
    const i = wildType().heavy.kabat.indexOf('33');
    const section = s.viewer.frame!.layout.sections[1];
    const x = (section.axis!.x(i) + section.axis!.x(i + 1)) / 2;
    const describer = () =>
      new Describer({
        spec: s.viewer.spec,
        store: s.viewer.store,
        addresses: s.viewer.addresses,
        view: s.viewer.view,
        layout: () => s.viewer.frame?.layout,
      });
    const cards: { title?: string; rows: readonly { label: string; value: string }[]; text: string }[] = [];
    for (const track of ['cdrs-kabat', 'cdrs-imgt']) {
      const t = s.viewer.frame!.layout.sections[1].tracks.find((k) => k.track.address.endsWith(`/${track}`))!;
      s.viewer.input.move(x, t.box.y + t.box.height / 2);
      await settle(s);
      const card = describer().describe(s.viewer.hover.current)!;
      expect(card.title).toBe('CDR-H1');
      cards.push(card);
    }
    expect(cards[0].rows).toContainEqual({ label: 'Track', value: 'CDRs (Kabat)' });
    expect(cards[1].rows).toContainEqual({ label: 'Track', value: 'CDRs (IMGT)' });
    expect(cards[0].text).not.toBe(cards[1].text);
  });
});
