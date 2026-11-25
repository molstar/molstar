// @vitest-environment jsdom
/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The spec's tooltips on the card (_spec/README.md "Hover card", Tooltips; the `seq:tooltips` provider): a row after
// the place's rows on the marks a note covers, none on an empty band or a gap inside it (EC-254) or from a hidden
// layer's tooltip (EC-255); an item said whole keeps its note (EC-260); a grid cell in range by x shows the row once,
// and a cell or a pair whose two parts meet one note counts it once (EC-258); several matches share the cap and the
// rest is counted, a density bin answers over its stretch (EC-259); values escaped and links Seq*'s own alone (EC-256);
// a failed table's issue (EC-174); the live region's words.

import { createSeqSpec } from '@molstar/seqx-core/lang';
import { DomHost } from '@molstar/ui/dom';
import { afterEach, beforeEach, describe, expect, it, vi } from 'vitest';
import { viewerItem } from '../config.js';
import { DefaultSeqViewerSpec } from '../default-spec.js';
import { describerFor } from '../describe.js';
import { cutItems } from '../hover-card.js';
import { createSeqViewer, type SeqViewer } from '../viewer.js';
import { type Frames, installFrames } from './frames.js';

let frames: Frames;
const viewers: SeqViewer[] = [];

beforeEach(() => {
  frames = installFrames();
  vi.spyOn(DomHost, 'now').mockImplementation(() => 1000);
});

afterEach(() => {
  for (const v of viewers.splice(0)) v.retire();
  frames.restore();
  vi.restoreAllMocks();
});

const Letters = 'MEEPQSDPSVEPPLSQETFSDLWKLLPENNVLSPLPSQAMDDLMLSPDDIEQWFTEDPGP';

type Build = ReturnType<typeof createSeqSpec>;

/** p53's start with domains (one hidden copy), a grid of it against itself, and the tooltips `tips` adds. */
function spec(tips: (b: Build, p: string, domains: string) => void) {
  const b = createSeqSpec({ title: 'card tooltips' });
  const p = b.data.sequence('p', { alphabet: 'protein', letters: Letters });
  const domains = b.data.ranges('domains', {
    target: p,
    rows: { start: [2, 40], end: [12, 50], id: ['a', 'b'], label: ['TAD', 'PRD'], hla: ['A*02:01', 'B*07:02'] },
  });
  const g = b.data.grid('g', { x: p, y: p, values: Float32Array.from({ length: Letters.length ** 2 }, (_, i) => i) });
  // a value at every position but 20 (a gap: no mark there)
  const values = b.data.values('vals', {
    target: p,
    values: Float32Array.from({ length: Letters.length }, (_, i) => (i === 20 ? Number.NaN : i)),
  });
  b.views.view('v', (v) => {
    v.section('s', { of: p }, (s) => {
      s.track('letters', { label: 'p' }).layer('l', { type: 'letters' });
      s.track('domains', { label: 'Domains' }).layer('d', { type: 'ranges', data: domains });
      s.track('hidden', { label: 'Hidden' }).layer('h', { type: 'ranges', data: domains, hidden: true });
      s.track('values', { label: 'Values' }).layer('vals', { type: 'values', data: values });
    });
    v.section('plot', { layout: 'grid', of: p, of_y: p, follow_y: 'none' })
      .track('plot')
      .layer('g', { type: 'grid', data: g });
  });
  tips(b, p, domains);
  return b.build();
}

async function make(tips: (b: Build, p: string, domains: string) => void, cardItems = 8) {
  const v = createSeqViewer({
    ...DefaultSeqViewerSpec(),
    ui: undefined,
    config: [[viewerItem('hover.cardItems'), cardItems]] as never,
  });
  viewers.push(v);
  await v.load(spec(tips));
  v.resize(900);
  await v.settle();
  frames.flush();
  v.draw();
  return v;
}

const P = '/data/p';

/** The card at a place of the main section: on a layer (by its track's index), an item, a position. */
function card(v: SeqViewer, track: number, position: number, element?: string) {
  const d = describerFor({ ...v, layout: () => v.frame?.layout });
  const s = v.frame!.layout.sections[0];
  const l = s.tracks[track].layers[0];
  return d.describe({
    section: s.section.address,
    layer: l.key,
    system: P,
    position,
    ...(element ? { element: { owner: '/data/domains', key: element } } : {}),
  });
}

const note = (b: Build, p: string) =>
  b.views.tooltip('notes', {
    target: p,
    rows: { start: [3], end: [6], tooltip: ['R175H'] },
    text: '{}, a **hotspot**',
    label: 'Note',
  });

describe('the spec’s tooltips on the card', () => {
  it('a row after the place’s rows on the letters a note covers, and in the live region’s words', async () => {
    const v = await make(note);
    const c = card(v, 0, 4)!;
    expect(c.rows.at(-1)).toEqual({
      label: 'Note',
      value: 'R175H, a hotspot',
      rich: [
        { kind: 'text', text: 'R175H, a ' },
        { kind: 'strong', children: [{ kind: 'text', text: 'hotspot' }] },
      ],
    });
    expect(c.rows.map((r) => r.label)).toEqual(['index', 'Letter', 'Values', 'TAD', 'Note']);
    expect(c.text).toContain('Note R175H, a hotspot');
    // outside it, no row
    expect(card(v, 0, 20)!.rows.some((r) => r.label === 'Note')).toBe(false);
  });

  it('none on an empty band or a gap inside the note (EC-254), none from a hidden layer’s tooltip (EC-255)', async () => {
    const wide = (b: Build, p: string) => {
      // a note over 14 to 30: the domains' band is empty there (TAD 2–12, PRD 40–50), the values' 20 is a gap
      b.views.tooltip('wide', { target: p, rows: { start: [14], end: [30], tooltip: ['wide'] }, label: 'Wide' });
      // by data on the hidden layer: its marks are not drawn, so no card says it
      b.views.tooltip('hidden', { on: '/views/v/s/hidden/h', text: 'from the hidden layer', label: 'Hidden' });
    };
    const v = await make(wide);
    const has = (c: ReturnType<typeof card>, label: string) => c?.rows.some((r) => r.label === label) ?? false;
    const values = v.frame!.layout.sections[0].tracks.findIndex((t) => t.track.label === 'Values');
    // inside the note: the letters answer, the values at 19 answer, the empty band and the value gap do not
    expect(has(card(v, 0, 20), 'Wide')).toBe(true);
    expect(has(card(v, values, 19), 'Wide')).toBe(true);
    expect(has(card(v, 1, 20), 'Wide')).toBe(false);
    expect(has(card(v, values, 20), 'Wide')).toBe(false);
    // the hidden layer draws nothing; the same items on the shown layer never carry its tooltip
    const hidden = v.frame!.layout.sections[0].tracks.find((t) => t.track.label === 'Hidden');
    expect(hidden?.layers.length ?? 0).toBe(0);
    expect(has(card(v, 1, 4, 'a'), 'Hidden')).toBe(false);
  });

  it('an item said whole keeps its note (EC-260)', async () => {
    const v = await make(note);
    const c = card(v, 1, 4, 'a')!;
    expect(c.title).toBe('TAD');
    expect(c.rows.at(-1)?.label).toBe('Note');
  });

  it('a grid cell in a note’s range by x alone shows its row once (EC-258)', async () => {
    const v = await make(note);
    const d = describerFor({ ...v, layout: () => v.frame?.layout });
    const s = v.frame!.layout.sections[1];
    const l = s.tracks[0]?.layers[0] ?? s.lays.flat()[0] ?? s.sectionLays[0];
    const c = d.describe({
      section: s.section.address,
      layer: l.key,
      system: P,
      position: 4,
      ySystem: P,
      yPosition: 30,
    });
    expect(c?.rows.filter((r) => r.label === 'Note')).toHaveLength(1);
    const both = d.describe({
      section: s.section.address,
      layer: l.key,
      system: P,
      position: 4,
      ySystem: P,
      yPosition: 5,
    });
    // both of the cell's positions in the note: one row, counted once (no `+1 more`, no `and 1 more`)
    expect(both?.rows.filter((r) => r.label === 'Note')).toHaveLength(1);
    expect(both?.moreCount).toBeUndefined();
    expect(both?.text).not.toContain('more');
  });

  it('several matches share the cap with the place’s items; the rest counted, `+k more` (EC-259)', async () => {
    const n = 30;
    const v = await make(
      (b, p) =>
        b.views.tooltip('many', {
          target: p,
          rows: {
            start: Array(n).fill(0),
            end: Array(n).fill(10),
            tooltip: Array.from({ length: n }, (_, i) => `n${i}`),
          },
        }),
      3,
    );
    const c = card(v, 0, 4)!;
    // the TAD item fills one of three: two tooltip rows made, the rest counted
    const tips = c.rows.filter((r) => r.label.startsWith('n'));
    expect(tips.map((r) => r.label)).toEqual(['n0', 'n1']);
    expect(c.moreCount).toBe(n - 2);
    expect(c.text).toContain(`and ${n - 2} more`);
    const drawn = cutItems(c, 3);
    expect(drawn.rows.at(-1)?.label).toBe(`+${n - 2} more`);
  });

  it('by data: values escaped, links Seq*’s own alone (EC-256)', async () => {
    const v = await make((b, _, domains) =>
      b.views.tooltip('hla', {
        on: domains,
        text: '{hla} [focus](!seq-focus=/refs/x) [out](javascript:alert(1)) [rel](/refs/dbd)',
        label: 'HLA',
      }),
    );
    const c = card(v, 1, 42, 'b')!;
    const row = c.rows.find((r) => r.label === 'HLA')!;
    expect(row.value).toBe('B*07:02 focus out rel');
    const links = JSON.stringify(row.rich);
    expect(links).toContain('!seq-focus=/refs/x');
    expect(links).not.toContain('javascript');
    expect(links).not.toContain('"href":"/refs/dbd"');
  });
});

describe('the spec’s tooltips on a density bin (EC-259)', () => {
  it('a zoomed-out bin answers over its stretch, not the pointer’s one position', async () => {
    const length = 200_000;
    const b = createSeqSpec({ title: 'bins' });
    const p = b.data.sequence('p', { alphabet: 'dna', letters: 'ACGT'.repeat(length / 4) });
    const n = 5000;
    const many = b.data.ranges('many', {
      target: p,
      rows: {
        start: Int32Array.from({ length: n }, (_, i) => i * 40),
        end: Int32Array.from({ length: n }, (_, i) => i * 40 + 20),
      },
    });
    // a note a few positions past the bin's first: the pointer's position alone misses it
    b.views.tooltip('note', { target: p, rows: { start: [30], end: [35], tooltip: ['in the bin'] }, label: 'Note' });
    b.views.view('v', (v) =>
      v.section('s', { of: p }, (s) => s.track('many', { label: 'Many' }).layer('m', { type: 'ranges', data: many })),
    );
    const v = createSeqViewer({ ...DefaultSeqViewerSpec(), ui: undefined });
    viewers.push(v);
    await v.load(b.build());
    v.resize(900);
    await v.settle();
    frames.flush();
    v.draw();
    const s = v.frame!.layout.sections[0];
    const l = s.tracks[0].layers[0];
    expect((l.lay as { mode?: string }).mode).toBe('density');
    const d = describerFor({ ...v, layout: () => v.frame?.layout });
    const c = d.describe({ section: s.section.address, layer: l.key, system: '/data/p', position: 0 });
    expect(c?.rows.find((r) => r.label === 'Note')?.value).toBe('in the bin');
  });
});

describe('a tooltip table failed under the mark (EC-174)', () => {
  it('gives no row and its issue on the tooltip, said once', async () => {
    const v = await make(note);
    const d = describerFor({ ...v, layout: () => v.frame?.layout });
    const said: string[] = [];
    d.compose(
      DefaultSeqViewerSpec().registry!.filter((e) => e.kind === 'card-provider') as never,
      8,
      () => {},
      (i) => said.push(`${i.code} ${i.address}`),
    );
    const s = v.frame!.layout.sections[0];
    const store = v.store.store;
    // every tile of the note's table reads as failed
    const get = store.get.bind(store);
    vi.spyOn(store, 'get').mockImplementation((k) => {
      const e = get(k);
      return e && k.includes('/views/notes') ? ({ ...e, status: 'failed' } as never) : e;
    });
    const place = { section: s.section.address, layer: s.tracks[0].layers[0].key, system: P, position: 4 };
    const c = d.describe(place);
    d.describe(place);
    expect(c?.rows.some((r) => r.label === 'Note')).toBe(false);
    expect(said).toEqual(['seq.viewer.card.tooltip /views/notes']);
  });
});

describe('an export’s card titles', () => {
  it('compose each item’s card from the export’s own layout, never the viewer’s frame', async () => {
    const v = await make(note);
    // the frame gone (a viewer between draws): the export's titles still say each item's card and its note
    Object.defineProperty(v, 'frame', { get: () => undefined, configurable: true });
    const out = await v.renderSvg({ idPrefix: 't-', titles: 'card' });
    expect(out.svg).toMatch(/<title>TAD[^<]*Note: R175H, a hotspot/);
  });
});
