// @vitest-environment jsdom
/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The hover card (_spec/README.md "Behaviors", hoverCard; "Ports", Hover card) and the keys panel: the card's content
// is the announcer's description of the place (a letter in its numberings, an item, a value, a grid cell, a pair, a
// row, a block's copy), headless on the HoverModel; its place is beside the pointer, never over the hovered box and
// inside the view; it shows on rest and hides on a drag, on Esc and on leave; a rectangle dragged on a grid shows its
// summary; `hover.card` off removes it. The keys panel lists every bound command once with its keys.

import { createSeqSpec } from '@molstar/seqx-core/lang';
import type { Hit as RenderHit } from '@molstar/seqx-core/render';
import { DomHost } from '@molstar/ui/dom';
import { formatChord } from '@molstar/ui/commands/keys';
import { afterEach, beforeEach, describe, expect, it, vi } from 'vitest';
import { Announcer } from '../announcer.js';
import { parseChord } from '../bindings.js';
import { RecordingContext } from '../canvas/_test/recording.js';
import { viewerItem } from '../config.js';
import { DefaultSeqViewerSpec } from '../default-spec.js';
import { Describer } from '../describe.js';
import { type CardBox, type HoverCard, type HoverCardPlace, placeCard } from '../hover-card.js';
import { describeKeys } from '../keys.js';
import { HoverCardProvider } from '../behaviors/hover-card.js';
import { KeysPanelProvider } from '../behaviors/keys-panel.js';
import { createSeqViewer, type SeqViewer } from '../viewer.js';
import { type Frames, installFrames } from './frames.js';

let frames: Frames;
const viewers: SeqViewer[] = [];
let clock = 1000;

beforeEach(() => {
  frames = installFrames();
  clock = 1000;
  vi.spyOn(DomHost, 'now').mockImplementation(() => clock);
  mockCanvas();
});

afterEach(() => {
  for (const v of viewers.splice(0)) v.retire();
  frames.restore();
  vi.restoreAllMocks();
  document.body.innerHTML = '';
});

const H = 'EVQLVESGGGLVQPGGSLRL';
// biome-ignore format: a label a position, in index order
const Kabat = [null, '1', '2', '3', '4', '5', '52', '52A', '53', '54', '82', '82A', '82B', '82C', '83', '100', '100A', '100B', '101', '102'];
// biome-ignore format: a label a position, in index order
const Imgt = ['1', '2', '3', '4', '5', '6', '55', '56', '57', '58', '105', '106', '107', '108', '109', '111', '111.1', '112.1', '112', '113'];

/** A heavy chain in two numberings with its domains, a value track and a pair; a grid of it against itself. */
function cardSpec(domain = 'CDR1') {
  const b = createSeqSpec({ title: 'card' });
  const h = b.data.sequence('h', { alphabet: 'protein', letters: H }, (s) => {
    s.numbering('kabat', { method: 'labels', labels: Kabat, scheme: 'kabat' });
    s.numbering('imgt', { method: 'labels', labels: Imgt, scheme: 'imgt' });
  });
  const d = b.data.ranges('d', {
    target: h,
    rows: { start: [2], end: [8], id: ['A'], label: [domain], category: ['loop'] },
  });
  const v = b.data.values('v', { target: h, values: Float32Array.from({ length: 20 }, (_, i) => i / 4) });
  const p = b.data.pairs('p', { a: h, rows: { a_start: [2], a_end: [3], b_start: [15], b_end: [16], id: ['2-15'] } });
  const g = b.data.grid('g', { x: h, y: h, values: Float32Array.from({ length: 400 }, (_, i) => i / 10) });
  b.views.view('v', (vw) => {
    vw.section('s', { of: h }, (s) => {
      s.track('ruler', { label: 'kabat' }).layer('r', { type: 'ruler', numbering: 'kabat' });
      s.track('letters', { label: 'h' }).layer('l', { type: 'letters' });
      s.track('domains', { label: 'Domains' }).layer('d', { type: 'ranges', data: d });
      s.track('values', { label: 'score' }).layer('v', { type: 'values', data: v });
      s.track('pairs', { label: 'contacts' }).layer('p', { type: 'pairs', data: p });
    });
    vw.section('plot', { layout: 'grid', of: h, of_y: h, follow_y: 'none' })
      .track('plot')
      .layer('g', { type: 'grid', data: g });
  });
  return b.build();
}

/** An alignment's rows. */
function rowsSpec() {
  const b = createSeqSpec({ title: 'rows' });
  const aln = b.data.alignment('aln', {
    alphabet: 'protein',
    rows: { label: ['human', 'mouse', 'zebrafish'], letters: ['MEEPQ', 'MEESQ', 'MDEPQ'] },
  });
  b.views.view('v', (v) =>
    v.section('s', { of: aln }).track('rows', { height: 60 }).layer('r', { type: 'rows', data: aln }),
  );
  return b.build();
}

/** A complex of two CDK2 copies and a cyclin, the copies split: a block a copy. */
function blocksSpec() {
  const b = createSeqSpec({ title: 'blocks' });
  const a = b.data.sequence('cdk2', { alphabet: 'protein', letters: 'ACDEFGHIKLMNPQRSTVWY'.repeat(3) });
  const c = b.data.sequence('cyclin', { alphabet: 'protein', letters: 'YWVTSRQPNMLKIHGFEDCA'.repeat(2) });
  b.views.view('v', (v) => {
    const s = v.section('complex', { copies: 'each' });
    s.block('a', { of: a, label: 'CDK2' });
    s.block('b', { of: c, label: 'cyclin A' });
    s.block('c', { of: a, label: 'CDK2' });
    s.track('letters').layer('letters', { type: 'letters' });
  });
  return b.build();
}

/** A recording port: what the card was told, in order. */
function recorder() {
  const calls: ({ show: HoverCard; place: HoverCardPlace } | 'hide')[] = [];
  return {
    calls,
    last: () => calls.at(-1),
    port: () => ({
      show: (card: HoverCard, place: HoverCardPlace) => void calls.push({ show: card, place }),
      hide: () => void calls.push('hide'),
    }),
  };
}

async function make(
  spec: ReturnType<typeof cardSpec>,
  options: { port?: ReturnType<typeof recorder>; config?: [unknown, unknown][]; target?: HTMLElement } = {},
) {
  const base = DefaultSeqViewerSpec();
  const v = createSeqViewer(
    {
      ...base,
      // mounted: the canvas alone (its parts are the behaviors'), else headless
      ui: options.target ? { canvas: true } : undefined,
      behaviors: base.behaviors.map((b) =>
        b.provider === HoverCardProvider && options.port
          ? { provider: HoverCardProvider, params: { port: options.port.port } }
          : b,
      ),
      config: [[viewerItem('viewport.maxPxPerPosition'), 1000], ...(options.config ?? [])] as never,
    },
    options.target,
  );
  viewers.push(v);
  await v.load(spec);
  v.resize(900);
  await v.settle();
  frames.flush();
  v.draw();
  return v;
}

const words = (v: SeqViewer) =>
  new Describer({ spec: v.spec, store: v.store, addresses: v.addresses, view: v.view, layout: () => v.frame?.layout });

/** The layer of a type in a section of the last frame, and the section. */
function layerOf(v: SeqViewer, section: number, type: string) {
  const s = v.frame!.layout.sections[section];
  for (const t of s.tracks) for (const l of t.layers) if (l.type.name === type) return { s, t, l };
  for (const ls of s.lays) for (const l of ls) if (l.type.name === type) return { s, t: undefined, l };
  for (const l of s.sectionLays) if (l.type.name === type) return { s, t: undefined, l };
  throw new Error(`no ${type} layer`);
}

/** The card for a hit set on the hover model, and the announcer's words for it. */
function cardOf(v: SeqViewer, hit: Partial<RenderHit> & { system: string }) {
  const owner = {};
  v.hover.set(hit as never, owner);
  const card = words(v).describe(v.hover.current);
  const a = new Announcer(v);
  const said = a.hoverText();
  a.dispose();
  v.hover.clear(owner);
  return { card, said };
}

const Heavy = '/data/h';

/** Canvases draw into recording contexts (jsdom has none). */
function mockCanvas() {
  vi.spyOn(HTMLCanvasElement.prototype, 'getContext').mockImplementation(function (this: HTMLCanvasElement) {
    return ((this as never as { __ctx?: RecordingContext }).__ctx ??= new RecordingContext()) as never;
  });
}

describe('the hover card’s content, headless on the HoverModel', () => {
  it('names a letter in every numbering of its system, and says what the announcer says', async () => {
    const v = await make(cardSpec());
    const { s, l } = layerOf(v, 0, 'letters');
    const { card, said } = cardOf(v, { layer: l.key, section: s.section.address, system: Heavy, position: 12 });
    expect(card).toEqual({
      title: 'h 13',
      rows: [
        { label: 'index', value: '12' },
        { label: 'kabat', value: '82B' },
        { label: 'imgt', value: '107' },
        { label: 'Letter', value: 'Q' },
        { label: 'score', value: '3' },
      ],
      text: 's, h 13 (kabat 82B, imgt 107), letter Q, score 3',
    });
    expect(said).toBe(card!.text);
  });

  it('names an item across a cyclic origin by its two runs from its start, the card and the live region alike', async () => {
    const b = createSeqSpec({ title: 'ring' });
    const ring = b.data.sequence('ring', { alphabet: 'dna', letters: 'ACGT'.repeat(25), cyclic: true, name: 'pR' });
    const f = b.data.ranges('f', {
      target: ring,
      rows: { start: [90], end: [10], id: ['bla'], label: ['bla (AmpR)'] },
    });
    b.views.view('main', (view) => view.section('s', { of: ring }).track('t').layer('l', { type: 'ranges', data: f }));
    const v = await make(b.build() as never);
    const item = layerOf(v, 0, 'ranges');
    const { card, said } = cardOf(v, {
      layer: item.l.key,
      section: item.s.section.address,
      system: '/data/ring',
      element: { owner: item.l.layer.data!.address, key: 'bla' },
    });
    expect(card?.rows).toContainEqual({ label: 'Extent', value: 'pR 91–100, 1–10' });
    expect(said).toBe(card!.text);
    expect(said).toContain('bla (AmpR), pR 91–100, 1–10');
  });

  it('a place’s items are those covering it: an item across the origin is not listed far from both its runs', async () => {
    const b = createSeqSpec({ title: 'ring' });
    const ring = b.data.sequence('ring', { alphabet: 'dna', letters: 'ACGT'.repeat(25), cyclic: true, name: 'pR' });
    const f = b.data.ranges('f', {
      target: ring,
      rows: { start: [90, 40], end: [10, 60], id: ['bla', 'ori'], label: ['bla (AmpR)', 'ori'] },
    });
    b.views.view('main', (view) => view.section('s', { of: ring }).track('t').layer('l', { type: 'ranges', data: f }));
    const v = await make(b.build() as never);
    const { s, l } = layerOf(v, 0, 'ranges');
    const at = (position: number) =>
      cardOf(v, { layer: l.key, section: s.section.address, system: '/data/ring', position }).card!.rows.map(
        (r) => r.label,
      );
    // at 50: ori alone (bla runs 91–100 and 1–10)
    expect(at(50)).toContain('ori');
    expect(at(50)).not.toContain('bla (AmpR)');
    // on either run: bla
    expect(at(95)).toContain('bla (AmpR)');
    expect(at(5)).toContain('bla (AmpR)');
  });

  it('says not numbered by a numbering with no label at the place', async () => {
    const v = await make(cardSpec());
    const { s, l } = layerOf(v, 0, 'letters');
    // Kabat has no label at the first position
    const { card } = cardOf(v, { layer: l.key, section: s.section.address, system: Heavy, position: 0 });
    expect(card?.rows).toContainEqual({ label: 'kabat', value: 'not numbered' });
    expect(card?.text).toBe('s, h 1 (imgt 1), not numbered by kabat, letter E, score 0');
  });

  it('names an item by its label, extent and category; a value track’s value; a pair’s two ends', async () => {
    const v = await make(cardSpec());
    const item = layerOf(v, 0, 'ranges');
    const owner = item.l.layer.data!.address;
    const a = cardOf(v, {
      layer: item.l.key,
      section: item.s.section.address,
      system: Heavy,
      element: { owner, key: 'A' },
    });
    expect(a.card).toEqual({
      title: 'CDR1',
      rows: [
        { label: 'Extent', value: 'h 3–8' },
        { label: 'Category', value: 'loop' },
      ],
      text: 'CDR1, h 3–8, loop',
    });
    expect(a.said).toBe(a.card!.text);
    // inside the item, the value track names its value at the position
    const val = layerOf(v, 0, 'values');
    const b = cardOf(v, { layer: val.l.key, section: val.s.section.address, system: Heavy, position: 6 });
    expect(b.card!.rows).toContainEqual({ label: 'score', value: '1.5' });
    expect(b.card!.rows).toContainEqual({ label: 'CDR1', value: 'h 3–8', item: true });
    expect(b.said).toBe(b.card!.text);
    const pair = layerOf(v, 0, 'pairs');
    const c = cardOf(v, {
      layer: pair.l.key,
      section: pair.s.section.address,
      system: Heavy,
      element: { owner: pair.l.layer.data!.address, key: '2-15' },
    });
    expect(c.card).toEqual({
      title: '2-15',
      rows: [
        { label: 'From', value: 'h 3' },
        { label: 'To', value: 'h 16' },
      ],
      text: '2-15, h 3 to h 16',
    });
    expect(c.said).toBe(c.card!.text);
  });

  it('names a grid cell by its two residues (letter and position) and its value', async () => {
    const v = await make(cardSpec());
    const g = layerOf(v, 1, 'grid');
    const { card, said } = cardOf(v, {
      layer: g.l.key,
      section: g.s.section.address,
      system: Heavy,
      position: 3,
      ySystem: Heavy,
      yPosition: 5,
      target: 'cell',
    });
    expect(card).toEqual({
      title: 'h L4 × h E6',
      rows: [
        { label: 'x', value: 'h L4' },
        { label: 'y', value: 'h E6' },
        { label: 'Value', value: '6.5' },
      ],
      text: 'plot, h L4 × h E6, value 6.5',
    });
    expect(said).toBe(card!.text);
  });

  it('names a row by its label', async () => {
    const v = await make(rowsSpec());
    const r = layerOf(v, 0, 'rows');
    const rows = (r.l.lay as { rowsSystem: string }).rowsSystem;
    const { card, said } = cardOf(v, {
      layer: r.l.key,
      section: r.s.section.address,
      system: rows,
      position: 1,
      target: 'row',
    });
    // the title names the row once: no row repeats it
    expect(card?.title).toBe('mouse');
    expect(card?.rows.map((r) => r.value)).not.toContain('mouse');
    expect(card?.text.startsWith('row mouse')).toBe(true);
    expect(said).toBe(card!.text);
  });

  it('names a block’s chain and which copy it is', async () => {
    const v = await make(blocksSpec());
    const { s, l } = layerOf(v, 0, 'letters');
    const { card, said } = cardOf(v, {
      layer: l.key,
      section: s.section.address,
      system: '/data/cdk2',
      position: 4,
      block: 2,
    });
    expect(card?.rows).toContainEqual({ label: 'Chain', value: 'CDK2' });
    expect(card?.rows).toContainEqual({ label: 'Copy', value: 'copy 2 of 2' });
    expect(card?.text.endsWith('chain CDK2, copy 2 of 2')).toBe(true);
    expect(said).toBe(card!.text);
  });
});

/** A point on the letters of the first section at a position. */
function letterPoint(v: SeqViewer, position: number): [number, number] {
  const s = v.frame!.layout.sections[0];
  const t = s.tracks[1];
  return [(s.axis!.x(position) + s.axis!.x(position + 1)) / 2, t.box.y + t.box.height / 2];
}

const overlaps = (a: CardBox, b: CardBox) =>
  a.x < b.x + b.width && b.x < a.x + a.width && a.y < b.y + b.height && b.y < a.y + a.height;

describe('the hover card’s place and its life', () => {
  it('shows on rest beside the pointer, never over the hovered cell, inside the view', async () => {
    const r = recorder();
    const v = await make(cardSpec(), { port: r });
    const [x, y] = letterPoint(v, 12);
    v.input.move(x, y);
    // a move says nothing: the card waits for the rest
    expect(r.calls.filter((c) => c !== 'hide')).toHaveLength(0);
    clock += 500;
    frames.timers.fire();
    const shown = r.last();
    if (!shown || shown === 'hide') throw new Error('no card');
    expect(shown.show.title).toBe('h 13');
    const { anchor, view } = shown.place;
    // the anchor is the cell: it holds the pointer, the card is clear of it and in the view
    expect(x >= anchor.x && x <= anchor.x + anchor.width && y >= anchor.y && y <= anchor.y + anchor.height).toBe(true);
    const card = { x: shown.place.x, y: shown.place.y, width: 200, height: 96 };
    expect(overlaps(card, anchor)).toBe(false);
    expect(card.x >= view.x && card.y >= view.y && card.x + card.width <= view.x + view.width).toBe(true);
    expect(card.y + card.height <= view.y + view.height).toBe(true);
  });

  it('lists `hover.cardItems` items at most, then `+k more`; the shown state and `more` keep every one', async () => {
    const b = createSeqSpec({ title: 'many' });
    const seq = b.data.sequence('s', { alphabet: 'dna', letters: 'ACGT'.repeat(25) });
    const n = 12;
    const f = b.data.ranges('f', {
      target: seq,
      rows: {
        start: Array.from({ length: n }, (_, i) => i),
        end: Array.from({ length: n }, () => 60),
        id: Array.from({ length: n }, (_, i) => `f${i}`),
      },
    });
    b.views.view('main', (view) => {
      const sec = view.section('s', { of: seq });
      sec.track('letters').layer('letters', { type: 'letters' });
      sec.track('t').layer('l', { type: 'ranges', data: f });
    });
    const r = recorder();
    const v = await make(b.build() as never, { port: r, config: [[viewerItem('hover.cardItems'), 3]] });
    // the pointer on the letter at 30: the place's card, every feature under it
    const sec = v.frame!.layout.sections[0];
    const lt = sec.tracks[0].box;
    v.input.move((sec.axis!.x(30) + sec.axis!.x(31)) / 2, lt.y + lt.height / 2);
    clock += 500;
    frames.timers.fire();
    const shown = r.last();
    if (!shown || shown === 'hide') throw new Error('no card');
    const items = shown.show.rows.filter((x) => x.item);
    expect(items).toHaveLength(3);
    expect(shown.show.rows.at(-1)).toEqual({ label: `+${n - 3} more`, value: '' });
    expect(shown.show.more).toHaveLength(n - 3);
    // the state holds the whole card; the live region says every item
    const whole = v.behavior(HoverCardProvider)!.shown.value!.card;
    expect(whole.rows.filter((x) => x.item)).toHaveLength(n);
    expect(whole.text).toContain('f11');
  });

  it('placeCard: never over the anchor and inside the view, wherever there is room', () => {
    const view = { x: 0, y: 0, width: 600, height: 400 };
    const size = { width: 180, height: 80 };
    let seed = 7;
    const rand = () => ((seed = (seed * 16807) % 2147483647) - 1) / 2147483646;
    for (let i = 0; i < 200; i++) {
      const w = 4 + rand() * 120;
      const h = 4 + rand() * 120;
      const anchor = { x: rand() * (view.width - w), y: rand() * (view.height - h), width: w, height: h };
      const pointer = { x: anchor.x + rand() * w, y: anchor.y + rand() * h };
      const p = placeCard(anchor, pointer, view, size);
      const card = { x: p.x, y: p.y, ...size };
      expect(card.x).toBeGreaterThanOrEqual(view.x);
      expect(card.y).toBeGreaterThanOrEqual(view.y);
      expect(card.x + card.width).toBeLessThanOrEqual(view.x + view.width);
      expect(card.y + card.height).toBeLessThanOrEqual(view.y + view.height);
      expect(overlaps(card, anchor), JSON.stringify({ anchor, p })).toBe(false);
    }
    // above first, else below
    expect(placeCard({ x: 100, y: 200, width: 10, height: 10 }, { x: 105, y: 205 }, view, size).side).toBe('above');
    expect(placeCard({ x: 100, y: 20, width: 10, height: 10 }, { x: 105, y: 25 }, view, size).side).toBe('below');
  });

  it('hides while a drag pans, on Esc and on leave; the next rest shows it again', async () => {
    const r = recorder();
    const v = await make(cardSpec(), { port: r });
    const card = () => v.behavior(HoverCardProvider)!.shown.value;
    const [x, y] = letterPoint(v, 12);
    const rest = () => {
      clock += 500;
      frames.timers.fire();
    };
    v.input.move(x, y);
    rest();
    expect(card()).toBeDefined();
    // a pan drag on the ruler (an axis): hidden at the press, hidden while it moves, never shown during it
    const ruler = v.frame!.layout.sections[0].tracks[0].box;
    v.input.move(x, ruler.y + ruler.height / 2);
    rest();
    expect(card()).toBeDefined();
    v.input.down(x, ruler.y + ruler.height / 2);
    expect(card()).toBeUndefined();
    v.input.move(x + 20, ruler.y + ruler.height / 2);
    v.input.move(x + 40, ruler.y + ruler.height / 2);
    rest();
    expect(card()).toBeUndefined();
    v.input.up();
    v.input.move(x, y);
    rest();
    expect(card()).toBeDefined();
    // Esc hides it
    v.input.key({ key: 'Escape' });
    expect(card()).toBeUndefined();
    rest();
    expect(card()).toBeUndefined();
    // the pointer moves on and rests: shown; it leaves: hidden
    v.input.move(...letterPoint(v, 14));
    rest();
    expect(card()?.card.title).toBe('h 15');
    v.input.leave();
    expect(card()).toBeUndefined();
    expect(r.last()).toBe('hide');
  });

  it('a press within the rest clears it: no card follows', async () => {
    const r = recorder();
    const v = await make(cardSpec(), { port: r });
    const [x, y] = letterPoint(v, 12);
    v.input.move(x, y);
    clock += 200;
    v.input.down(x, y);
    v.input.up();
    clock += 1000;
    frames.timers.fire();
    expect(v.behavior(HoverCardProvider)!.shown.value).toBeUndefined();
    expect(r.calls.filter((c) => c !== 'hide')).toHaveLength(0);
  });

  it('no stale card: a zoom key re-places the cursor’s card, an edit under a resting card sends it away', async () => {
    const r = recorder();
    const v = await make(cardSpec(), { port: r });
    const card = () => v.behavior(HoverCardProvider)!.shown.value;
    // the keyboard alone (the pointer off the canvas: a hover of its own would replace the cursor's card)
    v.input.key({ key: 'ArrowRight' });
    v.input.key({ key: 'ArrowRight' });
    frames.flush();
    const before = card()!;
    const c = v.selection.cursor.value!;
    // a zoom: the frame is laid out again, the cell moved; the card follows it there
    expect(v.input.key({ key: '=' })).toBe(true);
    // the zoom is a jump: animated on the page's frames, timed by the clock
    for (let i = 0; i < 4; i++) {
      clock += 100;
      frames.flush();
    }
    const after = card()!;
    expect(after.card.title).toBe(`h ${c.position + 1}`);
    const s = v.frame!.layout.sections[0];
    const cell = Math.min(s.axis!.x(c.position), s.axis!.x(c.position + 1));
    expect(after.place.anchor.x).toBeCloseTo(cell, 5);
    expect(after.place.anchor.width).not.toBeCloseTo(before.place.anchor.width, 1);
    // an edit under a resting card: gone at the frame, back with the new picture's words after the rest
    const [x, y] = letterPoint(v, 12);
    v.input.move(x, y);
    clock += 500;
    frames.timers.fire();
    expect(card()?.card.title).toBe('h 13');
    await v.execute('layer.hide', { address: v.frame!.layout.sections[0].tracks[3].layers[0].key });
    frames.flush();
    expect(card()).toBeUndefined();
    clock += 500;
    frames.timers.fire();
    // the value track is hidden: the card no longer names its value
    expect(card()?.card.rows.some((row) => row.label === 'score')).toBe(false);
  });

  it('a zoom that re-hits the same item under a resting pointer keeps the cursor’s card (placed again)', async () => {
    const r = recorder();
    const v = await make(cardSpec(), { port: r });
    const card = () => v.behavior(HoverCardProvider)!.shown.value;
    // the pointer rests on the domain (positions 2 to 8), the keyboard cursor elsewhere shows its card
    const s = v.frame!.layout.sections[0];
    const t = s.tracks[2];
    v.input.move((s.axis!.x(5) + s.axis!.x(6)) / 2, t.box.y + t.box.height / 2);
    expect(v.hover.current.element?.key).toBe('A');
    v.input.key({ key: 'ArrowRight' });
    frames.flush();
    const c = v.selection.cursor.value!;
    expect(card()?.card.title).toBe(`h ${c.position + 1}`);
    const version = v.hover.version.value;
    // a zoom: the pointer's position moves within the domain, so the hover is hit again, the same item
    expect(v.input.key({ key: '=' })).toBe(true);
    for (let i = 0; i < 4; i++) {
      clock += 100;
      frames.flush();
    }
    expect(v.hover.version.value).not.toBe(version);
    expect(v.hover.current.element?.key).toBe('A');
    expect(card()?.card.title).toBe(`h ${c.position + 1}`);
  });

  it('a load that keeps the hovered item’s section, layer and key is a new hover: its card comes with the new words', async () => {
    const r = recorder();
    const v = await make(cardSpec(), { port: r });
    const card = () => v.behavior(HoverCardProvider)!.shown.value;
    const s = v.frame!.layout.sections[0];
    const t = s.tracks[2];
    v.input.move((s.axis!.x(5) + s.axis!.x(6)) / 2, t.box.y + t.box.height / 2);
    clock += 500;
    frames.timers.fire();
    expect(card()?.card.title).toMatch(/CDR1/);
    // Esc: the card goes until the hover changes
    v.input.key({ key: 'Escape' });
    expect(card()).toBeUndefined();
    // the data loaded again, the domain renamed: the pointer rests on the same item (key A) of another model
    await v.load(cardSpec('CDR-H1'));
    await v.settle();
    frames.flush();
    expect(v.hover.current.element?.key).toBe('A');
    clock += 500;
    frames.timers.fire();
    expect(card()?.card.title).toMatch(/CDR-H1/);
  });

  it('shows at the keyboard cursor when a key moves it', async () => {
    const r = recorder();
    const v = await make(cardSpec(), { port: r });
    v.input.move(...letterPoint(v, 3));
    v.input.key({ key: 'ArrowRight' });
    v.input.key({ key: 'ArrowRight' });
    // made after the frame that drew the cursor (a reveal may have moved the view), never from the frame before it
    expect(v.behavior(HoverCardProvider)!.shown.value).toBeUndefined();
    frames.flush();
    const shown = v.behavior(HoverCardProvider)!.shown.value;
    const c = v.selection.cursor.value!;
    expect(shown?.card.title).toBe(`h ${c.position + 1}`);
    expect(overlaps({ x: shown!.place.x, y: shown!.place.y, width: 200, height: 96 }, shown!.place.anchor)).toBe(false);
  });

  it('shows a rectangle’s summary on a grid while it is dragged', async () => {
    const r = recorder();
    const v = await make(cardSpec(), { port: r });
    // a cell's middle on the plot's two axes (its strip's)
    const strips = (
      v.placementOf('/views/v/plot') as unknown as {
        strips: { axes: { xOfUnit(u: number): number }[]; ys: { xOfUnit(u: number): number }[] };
      }
    ).strips;
    const at = (i: number, j: number): [number, number] => [
      strips.axes[0].xOfUnit(i + 0.5),
      strips.ys[0].xOfUnit(j + 0.5),
    ];
    v.input.move(...at(2, 2));
    v.input.down(...at(2, 2));
    v.input.move(...at(3, 3));
    v.input.move(...at(5, 4));
    frames.timers.fire();
    const shown = v.behavior(HoverCardProvider)!.shown.value;
    expect(v.selection.value.value.location.rects.length).toBeGreaterThan(0);
    expect(shown?.card.title).toBe('12 cells');
    expect(shown?.card.rows.map((row) => row.label)).toEqual(['Mean', 'Min', 'Max']);
    expect(shown?.card.text).toMatch(/^Rectangle: 12 cells, mean /);
    v.input.up();
    expect(v.behavior(HoverCardProvider)!.shown.value).toBeUndefined();
  });

  it('makes its port once and leaves a host’s own (disposeWithViewer false) to the host', async () => {
    const made: (HTMLElement | undefined)[] = [];
    let disposed = 0;
    const port = (container: HTMLElement | undefined) => {
      made.push(container);
      return { show() {}, hide() {}, disposeWithViewer: false, dispose: () => void disposed++ };
    };
    const base = DefaultSeqViewerSpec();
    const v = createSeqViewer({
      ...base,
      ui: undefined,
      behaviors: base.behaviors.map((b) =>
        b.provider === HoverCardProvider ? { provider: HoverCardProvider, params: { port } } : b,
      ),
    });
    viewers.push(v);
    await v.load(cardSpec());
    v.resize(900);
    await v.settle();
    frames.flush();
    v.input.move(...letterPoint(v, 12));
    clock += 500;
    frames.timers.fire();
    v.input.move(...letterPoint(v, 13));
    clock += 500;
    frames.timers.fire();
    // headless: made at the first card, once
    expect(made).toEqual([undefined]);
    // a mount after it: the first given up (not disposed: the host's), one made over the canvas
    const el = document.createElement('div');
    const detach = v.behavior(HoverCardProvider)!.attach(el);
    expect(made).toEqual([undefined, el]);
    detach();
    v.retire();
    expect(disposed).toBe(0);
  });

  it('is the default adapter’s card over the canvas, and `hover.card` off removes it', async () => {
    vi.stubGlobal('ResizeObserver', undefined);
    const target = document.createElement('div');
    document.body.append(target);
    const v = await make(cardSpec(), { target });
    v.input.move(...letterPoint(v, 12));
    clock += 500;
    frames.timers.fire();
    const el = target.querySelector<HTMLElement>('.seqx-hover-card')!;
    expect(el.classList.contains('msui-surface')).toBe(true);
    expect(el.getAttribute('aria-hidden')).toBe('true');
    expect(el.style.display).toBe('');
    expect(el.querySelector('.msui-header')!.textContent).toBe('h 13');
    expect([...el.querySelectorAll('.msui-label')].map((n) => n.textContent)).toContain('kabat');
    v.input.leave();
    expect(el.style.display).toBe('none');
    // off: no card, and its part over a canvas adds no element (headless, the mount's call made by hand)
    const off = document.createElement('div');
    const w = await make(cardSpec(), { config: [[viewerItem('hover.card'), false]] });
    const detach = w.behavior(HoverCardProvider)!.attach(off);
    w.input.move(...letterPoint(w, 12));
    clock += 500;
    frames.timers.fire();
    expect(off.querySelector('.seqx-hover-card')).toBeNull();
    expect(w.behavior(HoverCardProvider)!.shown.value).toBeUndefined();
    detach();
  });
});

describe('the keys panel', () => {
  it('lists every bound command once with its keys, grouped, and every gesture', async () => {
    const v = await make(cardSpec());
    const groups = describeKeys(v, false);
    expect(groups.map((g) => g.id)).toEqual(['move', 'zoom', 'select', 'edit', 'view']);
    const entries = groups.flatMap((g) => g.entries);
    for (const c of v.commands) {
      const listed = v.config.view.keys[c.id];
      const first = listed ? parseChord(listed[0] ?? '') : c.shortcut;
      const of = entries.filter((e) => e.action === c.id);
      if (!first) continue;
      expect(of, c.id).toHaveLength(1);
      expect(of[0].keys, c.id).toContain(formatChord(first, false));
    }
    expect(entries.find((e) => e.action === 'help.keys')?.keys).toEqual(['?']);
    expect(entries.find((e) => e.action === 'pan')?.gestures).toContain('Drag on an axis');
    expect(entries.find((e) => e.action === 'view.zoom-in')?.keys).toEqual(['=', '+']);
  });

  it('opens with ?, draws the list, and closes with Esc', async () => {
    vi.stubGlobal('ResizeObserver', undefined);
    const target = document.createElement('div');
    document.body.append(target);
    const v = await make(cardSpec(), { target });
    // its DOM a lazy chunk, asked for at the first open
    expect(target.querySelector('.seqx-keys-panel')).toBeNull();
    expect(v.input.key({ key: '?', shiftKey: true })).toBe(true);
    await v.settle();
    const panel = target.querySelector<HTMLElement>('.seqx-keys-panel')!;
    expect(panel.style.display).toBe('');
    expect(panel.getAttribute('aria-label')).toBe('Keys and gestures');
    expect(panel.textContent).toContain('Zoom in');
    // Esc closes it, and nothing else (a selection stays)
    v.addresses.select('/data/h/(range :from 1 :to 4)');
    expect(v.input.key({ key: 'Escape' })).toBe(true);
    expect(panel.style.display).toBe('none');
    expect(v.selection.value.value.location.parts.length).toBeGreaterThan(0);
  });
  it('keeps its keys: Enter on Close closes it, arrows move no cursor; the focus comes back to what opened it', async () => {
    vi.stubGlobal('ResizeObserver', undefined);
    const target = document.createElement('div');
    document.body.append(target);
    const v = await make(cardSpec(), { target });
    const container = target.querySelector<HTMLElement>('[role="application"]')!;
    const key = (el: Element, k: string) =>
      el.dispatchEvent(new KeyboardEvent('keydown', { key: k, bubbles: true, cancelable: true }));
    // opened from the canvas: its Close button takes the focus (once its DOM landed)
    container.focus();
    key(container, '?');
    await v.settle();
    const panel = target.querySelector<HTMLElement>('.seqx-keys-panel')!;
    expect(panel.style.display).toBe('');
    const close = panel.querySelector('button')!;
    expect(document.activeElement).toBe(close);
    // an arrow in the panel stays there: the cursor does not move
    key(close, 'ArrowRight');
    expect(v.selection.cursor.value).toBeUndefined();
    // Enter on Close closes it; the focus is the canvas's again
    key(close, 'Enter');
    expect(panel.style.display).toBe('none');
    expect(document.activeElement).toBe(container);
    // opened from another control (a toolbar's button): Esc on the panel closes it and gives that control the focus
    const button = document.createElement('button');
    document.body.append(button);
    button.focus();
    await v.execute('help.keys', {});
    // drawn again as it opens: its Close button is a new one
    const again = panel.querySelector('button')!;
    expect(document.activeElement).toBe(again);
    key(again, 'Escape');
    expect(panel.style.display).toBe('none');
    expect(document.activeElement).toBe(button);
  });

  it('Esc cancels a running drag before it closes the panel', async () => {
    const v = await make(cardSpec());
    const keys = v.behavior(KeysPanelProvider)!;
    await v.execute('help.keys', {});
    expect(keys.open.value).toBe(true);
    // a drag-select on the letters, running
    const [x, y] = letterPoint(v, 4);
    v.input.move(x, y);
    v.input.down(x, y);
    v.input.move(...letterPoint(v, 8));
    expect(v.input.dragging).toBe(true);
    expect(v.input.key({ key: 'Escape' })).toBe(true);
    expect(v.input.dragging).toBe(false);
    expect(keys.open.value).toBe(true);
    expect(v.input.key({ key: 'Escape' })).toBe(true);
    expect(keys.open.value).toBe(false);
  });
});
