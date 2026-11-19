/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// A cyclic sequence whose letters come from a source, in the viewer (_spec/README.md "SpecModel"): once the source is
// open the model is resolved again with the length its coverage gave, so a row wrapping the origin is an item and the
// ring is placed; a focus on an item across the origin takes its shortest arc, never the hull of its runs (EC-12).

import { inlineSeq, seqToValue } from '@molstar/seqx-core/data';
import { createSeqSpec } from '@molstar/seqx-core/lang';
import { afterEach, describe, expect, it } from 'vitest';
import { DefaultSeqViewerSpec } from '../default-spec.js';
import { createSeqViewer, type SeqViewer } from '../viewer.js';
import { type Frames, installFrames } from './frames.js';

const N = 600;
const Letters = 'ACGT'.repeat(N / 4);
const Ring = '/views/v/ring';

function spec(fromSource: boolean) {
  const b = createSeqSpec({ title: 'ring' });
  const file = b.data.source('file', {
    format: 'seq',
    value: seqToValue(
      inlineSeq(
        { name: 'c', length: N, alphabet: 'dna', cyclic: true },
        { letters: { form: 'letters', start: 0, letters: Letters } },
      ),
    ),
  });
  const c = b.data.sequence(
    'c',
    fromSource ? { alphabet: 'dna', source: file, cyclic: true } : { alphabet: 'dna', letters: Letters, cyclic: true },
  );
  const f = b.data.ranges('f', { target: c, rows: { start: [560, 100], end: [40, 200], id: ['across', 'inside'] } });
  b.views.view('v', (v) => {
    const ring = v.section('ring', { of: c, layout: 'rings' });
    ring.track('f').layer('f', { type: 'ranges', data: f });
    v.section('detail', { of: c, follow: ring }).track('l').layer('l', { type: 'letters' });
  });
  return b.build();
}

let frames: Frames | undefined;
const viewers: SeqViewer[] = [];

afterEach(() => {
  for (const v of viewers.splice(0)) v.retire();
  frames?.restore();
  frames = undefined;
});

async function make(fromSource: boolean) {
  frames = installFrames();
  const viewer = createSeqViewer({ ...DefaultSeqViewerSpec(), ui: undefined });
  viewers.push(viewer);
  await viewer.load(spec(fromSource));
  viewer.resize(600);
  await viewer.settle();
  frames.flush();
  return viewer;
}

describe('a cyclic sequence read from a source, in the viewer', () => {
  it('is resolved again with its source’s length: the wrapped row is an item, as over inline letters', async () => {
    const v = await make(true);
    const model = v.spec.model.value;
    expect(model.systems.get('/data/c')?.length).toBe(N);
    expect(model.issues).toEqual([]);
    expect(model.annotations.get('/data/f')?.items?.filter((i) => i.key === 'across').length).toBe(2);
    expect(v.frame?.layout.issues).toEqual([]);
  });

  it('a section of blocks over a source alone is placed: an attach that moves no entry still gives the lengths', async () => {
    frames = installFrames();
    const viewer = createSeqViewer({ ...DefaultSeqViewerSpec(), ui: undefined });
    viewers.push(viewer);
    const b = createSeqSpec({ title: 'blocks' });
    const file = b.data.source('file', {
      format: 'seq',
      value: seqToValue(
        inlineSeq(
          { name: 'c', length: N, alphabet: 'dna' },
          { letters: { form: 'letters', start: 0, letters: Letters } },
        ),
      ),
    });
    const c = b.data.sequence('c', { alphabet: 'dna', source: file });
    b.views.view('v', (v) => {
      const s = v.section('blocks', {});
      s.block('one', { of: c });
      s.block('two', { of: c, range: [0, 100] });
      s.track('l').layer('l', { type: 'letters' });
    });
    await viewer.load(b.build());
    viewer.resize(600);
    await viewer.settle();
    frames.flush();
    expect(viewer.view.resolved.value?.sections[0].axis.ranges).toEqual([0, N, 0, 100]);
    expect(viewer.frame?.layout.issues).toEqual([]);
  });

  for (const fromSource of [false, true])
    it(`a focus on the item across the origin takes its shortest arc (${fromSource ? 'a source' : 'inline'})`, async () => {
      const v = await make(fromSource);
      expect(v.addresses.focus('/data/f/(item :id across)')).toMatchObject({ done: true });
      const w = v.view.get(Ring)!;
      expect([w.start, w.end]).toEqual([560, 640]);
      // an item inside the ring is its own hull
      expect(v.addresses.focus('/data/f/(item :id inside)')).toMatchObject({ done: true });
      expect([v.view.get(Ring)!.start, v.view.get(Ring)!.end]).toEqual([100, 200]);
    });
});
