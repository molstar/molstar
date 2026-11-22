/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// A liftover across a bridge (core _spec "Render", A liftover; EC-20, EC-212): hover and selection carried through a
// map read by region both ways (an old gene lights its images on the new assembly, an image lights the old gene), and
// the card of an item with two images lists them ("2 places").

import { categoryColumn, type Node } from '@molstar/core/repr/tree';
import { binPairs, type PairsList, type SeqData, seqToValue } from '@molstar/seqx-core/data';
import { assembly, createSeqSpec } from '@molstar/seqx-core/lang';
import { type Location, location } from '@molstar/seqx-core/query';
import { type Hit, hitView } from '@molstar/seqx-core/render';
import { loadGenome } from '@molstar/seqx-core/resolve';
import { afterEach, beforeAll, describe, expect, it } from 'vitest';
import { DefaultSeqViewerSpec } from '../default-spec.js';
import { describerFor } from '../describe.js';
import { Projector } from '../models/projection.js';
import { createSeqViewer, type SeqViewer } from '../viewer.js';
import { type Frames, installFrames } from './frames.js';

let frames: Frames | undefined;
const viewers: SeqViewer[] = [];

beforeAll(async () => {
  await loadGenome();
});

afterEach(() => {
  for (const v of viewers.splice(0)) v.retire();
  frames?.restore();
  frames = undefined;
});

const Old = [{ name: 'chr1', length: 100_000 }];
const New = [
  { name: 'chr1', length: 80_000 },
  { name: 'chr2', length: 80_000 },
];
/** Old chr1's head straight, its 70–90 kb duplicated onto new chr1 and new chr2 (`a` old, `b` new). */
const Runs = [
  { from: 0, b: 'chr1', to: 0, length: 20_000 },
  { from: 70_000, b: 'chr1', to: 60_000, length: 20_000 },
  { from: 70_000, b: 'chr2', to: 60_000, length: 20_000 },
];

function source(): SeqData {
  const n = Runs.length;
  const list: PairsList = {
    form: 'list',
    aSystem: categoryColumn(Runs.map(() => 'chr1')),
    bSystem: categoryColumn(Runs.map((r) => r.b)),
    aStart: Int32Array.from(Runs, (r) => r.from),
    aEnd: Int32Array.from(Runs, (r) => r.from + r.length),
    bStart: Int32Array.from(Runs, (r) => r.to),
    bEnd: Int32Array.from(Runs, (r) => r.to + r.length),
    aStrand: new Int8Array(n).fill(1),
    bStrand: new Int8Array(n).fill(1),
    row: Int32Array.from(Runs, (_, i) => i),
  };
  const sides = [
    ['a', Old],
    ['b', New],
  ] as const;
  return {
    coverage: {
      systems: sides.flatMap(([side, cs]) => cs.map((c) => ({ ...c, side }))),
      shapes: ['chr1', 'chr2'].map((system) => ({ system, shape: 'pairs' as const })),
      levels: [1, 1024],
    },
    pairs: sides.flatMap(([side, cs]) =>
      cs.flatMap((c) =>
        [1, 1024].map((level) => ({
          request: { system: c.name, start: 0, end: c.length, level, shapes: ['pairs' as const], side },
          answer: { ok: true as const, tables: { pairs: level === 1 ? list : binPairs(list, level) } },
        })),
      ),
    ),
  };
}

const OldChr1 = '/data/old/(sequence :name chr1)';
const NewChr1 = '/data/new/(sequence :name chr1)';
const NewChr2 = '/data/new/(sequence :name chr2)';

function spec(): Node {
  const b = createSeqSpec({ title: 'a liftover' });
  const old = b.data.assembly('old', assembly('Old', Old, { order: 'source' } as never));
  const neu = b.data.assembly('new', assembly('New', New, { order: 'source' } as never));
  const src = b.data.source('chain', { format: 'seq', value: seqToValue(source()) as never });
  const map = b.data.map('lift', { from: old, to: neu, method: 'segments', source: src } as never);
  const genes = b.data.ranges('genes', {
    target: old,
    rows: { id: ['head', 'twice'], target: ['chr1', 'chr1'], start: [5_000, 72_000], end: [15_000, 78_000] },
  } as never);
  b.views.view('v', (v) => {
    v.section('old', { of: old } as never, (s) => s.track('genes').layer('genes', { type: 'ranges', data: genes }));
    v.bridge('br', { a: '/views/v/old', b: '/views/v/new' } as never, (br) =>
      br.layer('runs', { type: 'pairs', data: map } as never),
    );
    v.section('new', { of: neu } as never, (s) => s.track('genes').layer('genes', { type: 'ranges', data: genes }));
  });
  return b.build();
}

async function make(): Promise<SeqViewer> {
  frames = installFrames();
  const v = createSeqViewer(DefaultSeqViewerSpec());
  viewers.push(v);
  await v.load(spec());
  v.resize(900);
  await v.settle();
  frames.flush();
  v.draw();
  await v.settle();
  frames.flush();
  v.draw();
  return v;
}

/** The first hit on `section` naming the item `key` (a scan of the drawn page). */
function hitOn(v: SeqViewer, section: string, key: string, system: string): Hit {
  const layout = v.frame!.layout;
  const sec = layout.sections.find((s) => s.section.address === section)!;
  for (let y = sec.y; y < sec.y + sec.height; y += 2)
    for (let x = 0; x < 900; x += 2) {
      const h = hitView(layout, x, y);
      // as the input takes it: the section read from the layer
      if (h?.element?.key === key && h.system === system) return { ...h, section: undefined } as never as Hit;
    }
  throw new Error(`no hit on ${key}`);
}

const parts = (loc: Location | undefined) =>
  (loc?.parts ?? []).map((p) => [p.system, Array.from(p.ranges)] as const).sort((a, b) => a[0].localeCompare(b[0]));

describe('a liftover’s hover and selection, through the map both ways (EC-20, EC-212)', () => {
  it('an old gene lights its two images on the new assembly; its card lists them', async () => {
    const v = await make();
    const h = hitOn(v, '/views/v/old', 'twice', OldChr1);
    v.hover.set(h, 'test');
    const all = parts(v.hover.projectedAll());
    expect(all).toEqual([
      [NewChr1, [62_000, 68_000]],
      [NewChr2, [62_000, 68_000]],
      [OldChr1, [72_000, 78_000]],
    ]);
    const card = describerFor({
      spec: v.spec,
      store: v.store,
      addresses: v.addresses,
      view: v.view,
      layout: () => v.frame?.layout,
    }).describe(v.hover.current);
    expect(card?.rows).toContainEqual({ label: 'Places', value: '2 places: chr1 62,001–68,000, chr2 62,001–68,000' });
    expect(card?.text).toContain('in 2 places');
  });

  it('an image on the new assembly lights the old gene; one image has no places row', async () => {
    const v = await make();
    const h = hitOn(v, '/views/v/new', 'twice', NewChr2);
    v.hover.set(h, 'test');
    expect(parts(v.hover.projectedAll()).find((p) => p[0] === OldChr1)).toEqual([OldChr1, [72_000, 78_000]]);
    const card = describerFor({
      spec: v.spec,
      store: v.store,
      addresses: v.addresses,
      view: v.view,
      layout: () => v.frame?.layout,
    }).describe(v.hover.current);
    expect(card?.rows.some((r) => r.label === 'Places')).toBe(false);
  });

  it('a selection is carried both ways: the old gene to both images, an image back to the gene', async () => {
    const v = await make();
    const p = new Projector(v.spec, v.view, v.store.regions);
    const old = location({ parts: [{ system: OldChr1, ranges: [72_000, 78_000] }] });
    expect(parts(p.of(old))).toEqual([
      [NewChr1, [62_000, 68_000]],
      [NewChr2, [62_000, 68_000]],
    ]);
    const image = location({ parts: [{ system: NewChr1, ranges: [62_000, 63_000] }] });
    // one step: to the old gene (never on to the other image, which no map joins to new chr1)
    expect(parts(p.of(image))).toEqual([[OldChr1, [72_000, 73_000]]]);
  });
});
