// @vitest-environment jsdom
/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// A canonical transcript in the hover card and the announcer (_spec/README.md "Hover card"): a source's item whose
// `canonical` field says `Canonical` or `MANE Select` (GFF3's and GTF's tags) has a Transcript row saying so; an item
// with none says nothing.

import { categoryColumn } from '@molstar/core/repr/tree';
import type { FormatRuntime, RegionAnswer, RegionRequest } from '@molstar/seqx-core/data';
import { createSeqSpec } from '@molstar/seqx-core/lang';
import type { Hit as RenderHit } from '@molstar/seqx-core/render';
import { afterEach, beforeEach, describe, expect, it, vi } from 'vitest';
import { RecordingContext } from '../canvas/_test/recording.js';
import { Announcer } from '../announcer.js';
import { viewerItem } from '../config.js';
import { DefaultSeqViewerSpec } from '../default-spec.js';
import { Describer } from '../describe.js';
import { createSeqViewer, type SeqViewer } from '../viewer.js';
import { type Frames, installFrames } from './frames.js';

const N = 400;
const Letters = 'ACGT'.repeat(N / 4);

const Genes = [
  { id: 'mane', start: 20, end: 80, strand: -1, canonical: 'MANE Select' },
  { id: 'canon', start: 120, end: 180, strand: 1, canonical: 'Canonical' },
  { id: 'none', start: 220, end: 280, strand: 0, canonical: null },
];

/** A source read by region: the genes overlapping a region, with their strands. */
const Genome: FormatRuntime = {
  name: 'genome',
  formats: ['genome'],
  async open() {
    return {
      coverage: {
        systems: [{ name: 'chr', length: N, alphabet: 'dna' }],
        shapes: [
          { system: 'chr', shape: 'letters' },
          { system: 'chr', shape: 'ranges' },
        ],
        levels: [1],
      },
      read: async (r: RegionRequest): Promise<RegionAnswer> => {
        const tables: { letters?: unknown; ranges?: unknown } = {};
        if (r.shapes.includes('letters'))
          tables.letters = { form: 'letters', start: r.start, letters: Letters.slice(r.start, r.end) };
        if (r.shapes.includes('ranges')) {
          const gs = Genes.filter((g) => g.start < r.end && g.end > r.start);
          tables.ranges = {
            form: 'list',
            start: Int32Array.from(gs, (g) => g.start),
            end: Int32Array.from(gs, (g) => g.end),
            strand: Int8Array.from(gs, (g) => g.strand),
            id: categoryColumn(gs.map((g) => g.id)),
            fields: [{ name: 'canonical', column: categoryColumn(gs.map((g) => g.canonical)) }],
          };
        }
        return { ok: true, tables: tables as never };
      },
      dispose() {},
    };
  },
};

function spec() {
  const b = createSeqSpec({ title: 'genes' });
  const file = b.data.source('file', { format: 'genome', uri: 'genome:chr' });
  const c = b.data.sequence('chr', { alphabet: 'dna', source: file });
  const genes = b.data.ranges('genes', { target: c, source: file });
  b.views.view('v', (v) =>
    v.section('s', { of: c }).track('genes', { height: 40 }).layer('genes', { type: 'ranges', data: genes }),
  );
  return b.build();
}

let frames: Frames | undefined;
const viewers: SeqViewer[] = [];

beforeEach(() => {
  vi.spyOn(HTMLCanvasElement.prototype, 'getContext').mockImplementation(function (this: HTMLCanvasElement) {
    return ((this as never as { __ctx?: RecordingContext }).__ctx ??= new RecordingContext()) as never;
  });
});

afterEach(() => {
  vi.restoreAllMocks();
  for (const v of viewers.splice(0)) v.retire();
  frames?.restore();
  frames = undefined;
});

async function make() {
  frames = installFrames();
  const v = createSeqViewer({
    ...DefaultSeqViewerSpec(),
    ui: undefined,
    registry: [...(DefaultSeqViewerSpec().registry ?? []), { kind: 'format-runtime', runtime: Genome }],
    config: [[viewerItem('viewport.maxPxPerPosition'), 1000]] as never,
  });
  viewers.push(v);
  await v.load(spec());
  v.resize(600);
  await v.settle();
  frames.flush();
  for (let i = 0; i < 4; i++) {
    v.draw();
    await v.settle();
    frames.flush();
  }
  v.draw();
  return v;
}

function cardOf(v: SeqViewer, key: string) {
  const s = v.frame!.layout.sections[0];
  const l = s.tracks[0].layers[0];
  const owner = {};
  const hit: Partial<RenderHit> = {
    layer: l.key,
    section: s.section.address,
    system: l.layer.data!.system,
    element: { owner: l.layer.data!.address, key },
  };
  v.hover.set(hit as never, owner);
  const card = new Describer({
    spec: v.spec,
    store: v.store,
    addresses: v.addresses,
    view: v.view,
    layout: () => v.frame?.layout,
  }).describe(v.hover.current);
  const a = new Announcer(v);
  const said = a.hoverText();
  a.dispose();
  v.hover.clear(owner);
  return { card, said };
}

describe('a canonical transcript in the card and the announcement', () => {
  it('says Canonical or MANE Select as the Transcript row, as the data does; none says nothing', async () => {
    const v = await make();
    const mane = cardOf(v, 'mane');
    expect(mane.card!.rows).toContainEqual({ label: 'Transcript', value: 'MANE Select' });
    expect(mane.said).toBe(mane.card!.text);
    expect(mane.card!.text).toMatch(/, mane select$/);
    expect(cardOf(v, 'canon').card!.rows).toContainEqual({ label: 'Transcript', value: 'Canonical' });
    const none = cardOf(v, 'none');
    expect(none.card!.rows.map((r) => r.label)).not.toContain('Transcript');
    expect(none.said).toBe(none.card!.text);
  });
});
