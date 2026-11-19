// @vitest-environment jsdom
/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Declared numberings in the viewer (core _spec/README.md, Language, Numberings; _spec/README.md "Input",
// Announcements; "AddressModel"): a hover that rests says the position's label in the numbering the section's ruler
// shows, an address in a numbering is followed, a copied address is printed in the numbering asked; and the ruler's
// labels drawn on the canvas are the SVG's.

import { createSeqSpec } from '@molstar/seqx-core/lang';
import { rangesLocation } from '@molstar/seqx-core/query';
import { renderSvg } from '@molstar/seqx-core/render';
import { DomHost } from '@molstar/ui/dom';
import { afterEach, beforeEach, describe, expect, it, vi } from 'vitest';
import { Announcer } from '../announcer.js';
import { CanvasBackend } from '../canvas/backend.js';
import { FakeCanvas } from '../canvas/_test/recording.js';
import { viewerItem } from '../config.js';
import { DefaultSeqViewerSpec } from '../default-spec.js';
import { createSeqViewer, type SeqViewer } from '../viewer.js';
import { type Frames, installFrames } from './frames.js';

const H = 'EVQLVESGGGLVQPGGSLRL';
// biome-ignore format: a label a position, in index order
const Kabat = [null, '1', '2', '3', '4', '5', '52', '52A', '53', '54', '82', '82A', '82B', '82C', '83', '100', '100A', '100B', '101', '102'];
// biome-ignore format: a label a position, in index order
const Imgt = ['1', '2', '3', '4', '5', '6', '55', '56', '57', '58', '105', '106', '107', '108', '109', '111', '111.1', '112.1', '112', '113'];

let frames: Frames;
const viewers: SeqViewer[] = [];
let clock = 1000;

beforeEach(() => {
  frames = installFrames();
  clock = 1000;
  vi.spyOn(DomHost, 'now').mockImplementation(() => clock);
});

afterEach(() => {
  for (const v of viewers.splice(0)) v.retire();
  frames.restore();
  vi.restoreAllMocks();
  document.body.innerHTML = '';
});

/** A heavy chain shown twice: a section a numbering, its ruler over its letters. */
function spec() {
  const b = createSeqSpec({ title: 'numberings' });
  const h = b.data.sequence('h', { alphabet: 'protein', letters: H }, (s) => {
    s.numbering('kabat', { method: 'labels', labels: Kabat, scheme: 'kabat' });
    s.numbering('imgt', { method: 'labels', labels: Imgt, scheme: 'imgt' });
  });
  b.views.view('v', (v) => {
    for (const n of ['kabat', 'imgt'])
      v.section(n, { of: h }, (s) => {
        s.track('ruler', { label: n }).layer('r', { type: 'ruler', numbering: n });
        s.track('letters', { label: 'h' }).layer('l', { type: 'letters' });
      });
  });
  return b.build();
}

async function make(copyNumbering = 'ordinal') {
  const v = createSeqViewer({
    ...DefaultSeqViewerSpec(),
    ui: undefined,
    config: [
      [viewerItem('viewport.maxPxPerPosition'), 1000],
      [viewerItem('copyLink.numbering'), copyNumbering],
    ] as never,
  });
  viewers.push(v);
  await v.load(spec());
  v.resize(900);
  await v.settle();
  frames.flush();
  v.draw();
  return v;
}

describe('declared numberings in the viewer', () => {
  it('a hover that rests says the label in each numbering, the section’s first', async () => {
    const v = await make();
    const a = new Announcer(v, { minGapMs: 0, restMs: 500 });
    for (const [section, label] of [
      // every declared numbering, the section's shown one first
      [0, 'kabat 82B, imgt 107'],
      [1, 'imgt 107, kabat 82B'],
    ] as const) {
      const s = v.frame!.layout.sections[section];
      const t = s.tracks[1];
      const x = (s.axis!.x(12) + s.axis!.x(13)) / 2;
      v.input.move(x, t.box.y + t.box.height / 2);
      frames.flush();
      frames.timers.fire();
      expect(a.message.value.text).toBe(`${s.id}, h 13 (${label}), letter Q`);
      clock += 1000;
    }
    a.dispose();
  });

  it('an address in a numbering resolves; a copied one is printed in the numbering asked', async () => {
    const v = await make('kabat');
    const r = v.addresses.resolve('/data/h/(range :from 82 :to 83 :numbering kabat)');
    expect(r.status).toBe('resolved');
    expect(r.result?.type === 'location' && Array.from(r.result.location.parts[0].ranges)).toEqual([10, 15]);
    const loc = { type: 'location' as const, location: rangesLocation('/data/h', [11, 14]) };
    expect(v.addresses.addressOf(loc)).toBe('/data/h/(range :from 82A :numbering kabat :to 82C)');
  });
});

describe('a numbering from a source, with no ruler showing it', () => {
  it('a link naming it waits, asks for its labels, and resolves when they arrive', async () => {
    const b = createSeqSpec({ title: 'no ruler' });
    const t = b.data.source('kabat-table', {
      format: 'json',
      system: 'h',
      value: Kabat.map((label, position) => ({ position, label })),
    });
    const h = b.data.sequence('h', { alphabet: 'protein', letters: H }, (s) =>
      s.numbering('kabat', { method: 'labels', source: t, scheme: 'kabat' }),
    );
    b.views.view('v', (v) => v.section('s', { of: h }).track('l').layer('l', { type: 'letters' }));
    const v = createSeqViewer({ ...DefaultSeqViewerSpec(), ui: undefined });
    viewers.push(v);
    await v.load(b.build());
    v.resize(900);
    await v.settle();
    const address = '/data/h/(range :from 82 :to 83 :numbering kabat)';
    const first = v.addresses.resolve(address);
    expect(first.status).toBe('waiting');
    expect(first.issues.map((i) => i.code)).toEqual(['seq.query.numbering.pending', 'seq.address.waiting']);
    // the address model asked for the labels (no layer did)
    await Promise.resolve();
    await v.store.store.whenSettled();
    const r = v.addresses.resolve(address);
    expect(r.status).toBe('resolved');
    expect(r.result?.type === 'location' && Array.from(r.result.location.parts[0].ranges)).toEqual([10, 15]);
  });
});

describe('the ruler in a declared numbering: canvas and SVG', () => {
  it('draws the same labels', async () => {
    const r = await renderSvg(spec(), {
      size: { width: 1400, height: 240 },
      viewport: { start: 0, end: H.length },
      legend: 'none',
    });
    expect(r.issues).toEqual([]);
    const svg = [...r.svg.matchAll(/<text[^>]*>([^<]*)<\/text>/g)].map((m) => m[1]);
    const base = new FakeCanvas();
    const labelsCanvas = new FakeCanvas();
    const backend = new CanvasBackend(base, new FakeCanvas(), labelsCanvas, {
      createCanvas: () => new FakeCanvas(),
    });
    backend.resize(r.lists.width, r.lists.height, 1);
    backend.draw(r.lists, true);
    // the labels pass draws on the labels canvas
    const canvas = [...base.context.take(), ...labelsCanvas.context.take()]
      .filter((c) => c.startsWith('fillText('))
      .map((c) => JSON.parse(c.slice(9, c.indexOf('",') + 1)) as string);
    for (const l of ['52A', '82C', '100B', '111.1', '113']) {
      expect(svg).toContain(l);
      expect(canvas).toContain(l);
    }
    // every label of the ruler the SVG prints, the canvas draws
    const labels = new Set([...Kabat, ...Imgt].filter((l): l is string => l !== null));
    expect(canvas.filter((t) => labels.has(t)).sort()).toEqual(svg.filter((t) => labels.has(t)).sort());
  });
});
