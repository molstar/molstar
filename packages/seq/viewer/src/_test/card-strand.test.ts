// @vitest-environment jsdom
/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// A range item's strand in the hover card and the announcer (_spec/README.md "Hover card"): a gene on the minus strand
// of a source's ranges says `Strand −`, as an inline annotation's item does; one with no strand says none.

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
  { id: 'minus', start: 20, end: 80, strand: -1 },
  { id: 'plus', start: 120, end: 180, strand: 1 },
  { id: 'none', start: 220, end: 280, strand: 0 },
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
          };
        }
        return { ok: true, tables: tables as never };
      },
      dispose() {},
    };
  },
};

function spec(from: 'source' | 'inline') {
  const b = createSeqSpec({ title: 'genes' });
  const file = b.data.source('file', { format: 'genome', uri: 'genome:chr' });
  const c = b.data.sequence(
    'chr',
    from === 'source' ? { alphabet: 'dna', source: file } : { alphabet: 'dna', letters: Letters },
  );
  const genes =
    from === 'source'
      ? b.data.ranges('genes', { target: c, source: file })
      : b.data.ranges('genes', {
          target: c,
          rows: {
            start: Genes.map((g) => g.start),
            end: Genes.map((g) => g.end),
            id: Genes.map((g) => g.id),
            strand: ['-', '+', '.'],
          },
        });
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

async function make(from: 'source' | 'inline') {
  frames = installFrames();
  const v = createSeqViewer({
    ...DefaultSeqViewerSpec(),
    ui: undefined,
    registry: [...(DefaultSeqViewerSpec().registry ?? []), { kind: 'format-runtime', runtime: Genome }],
    config: [[viewerItem('viewport.maxPxPerPosition'), 1000]] as never,
  });
  viewers.push(v);
  await v.load(spec(from));
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

describe.each(['source', 'inline'] as const)('a range’s strand in the card and the announcement (%s)', (from) => {
  it('says a minus-strand gene’s strand, and a plus one’s; none when the data has none', async () => {
    const v = await make(from);
    const minus = cardOf(v, 'minus');
    expect(minus.card!.rows).toContainEqual({ label: 'Strand', value: '−' });
    expect(minus.card!.text).toBe('minus, chr 21–80, strand −');
    expect(minus.said).toBe(minus.card!.text);
    expect(cardOf(v, 'plus').card!.rows).toContainEqual({ label: 'Strand', value: '+' });
    const none = cardOf(v, 'none');
    expect(none.card!.rows.map((r) => r.label)).not.toContain('Strand');
    expect(none.said).toBe(none.card!.text);
  });
});
