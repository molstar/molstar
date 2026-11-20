/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The model's issues as the viewer says them (_spec/README.md "SpecModel", Issues of a pending address):
// `seq.resolve.length.unknown` of a row wrapping the origin of a cyclic sequence read from a source is held back while
// the source opens, and gone once its length arrives; one whose source failed is said; `seq.resolve.system.unnamed`
// (the source opened and names no system after the sequence) is said once it opened.

import type { Issue } from '@molstar/core/repr/tree';
import { inlineSeq, seqToValue } from '@molstar/seqx-core/data';
import { createSeqSpec } from '@molstar/seqx-core/lang';
import type { StaticModel } from '@molstar/seqx-core/resolve';
import { afterEach, describe, expect, it } from 'vitest';
import { DefaultSeqViewerSpec } from '../default-spec.js';
import { modelIssues, shownModelIssues } from '../model-issues.js';
import { createSeqViewer, type SeqViewer } from '../viewer.js';
import { type Frames, installFrames } from './frames.js';

const N = 600;
const Letters = 'ACGT'.repeat(N / 4);

/** A ring over a cyclic sequence read from a `seq` source holding system `c`, a row across its origin. */
function spec(id = 'c') {
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
  const c = b.data.sequence(id, { alphabet: 'dna', source: file, cyclic: true });
  const f = b.data.ranges('f', { target: c, rows: { start: [560], end: [40], id: ['across'] } });
  b.views.view('v', (v) =>
    v.section('ring', { of: c, layout: 'rings' }).track('f').layer('f', { type: 'ranges', data: f }),
  );
  return b.build();
}

let frames: Frames | undefined;
const viewers: SeqViewer[] = [];

afterEach(() => {
  for (const v of viewers.splice(0)) v.retire();
  frames?.restore();
  frames = undefined;
});

const codes = (issues: readonly Issue[]) => issues.map((i) => i.code);

describe('the model’s issues as the viewer says them', () => {
  it('`seq.resolve.length.unknown` is held back while its source opens, and gone once the length arrives', async () => {
    frames = installFrames();
    const v = createSeqViewer({ ...DefaultSeqViewerSpec(), ui: undefined });
    viewers.push(v);
    const loading = v.load(spec());
    // the model before the source opens: the row waits for the length, said by core, held back by the viewer
    expect(codes(v.spec.model.value.issues)).toContain('seq.resolve.length.unknown');
    expect(codes(modelIssues(v))).not.toContain('seq.resolve.length.unknown');
    await loading;
    v.resize(600);
    await v.settle();
    expect(codes(v.spec.model.value.issues)).toEqual([]);
    expect(modelIssues(v)).toEqual([]);
  });

  it('said when the source settled without the length (opened, or failed)', () => {
    const issue = { code: 'seq.resolve.length.unknown', severity: 'warning', message: '', address: '/data/f' };
    const model = {
      issues: [issue],
      maps: new Map(),
      annotations: new Map([['/data/f', { target: '/data/c' }]]),
      systems: new Map([['/data/c', { source: '/data/file' }]]),
    } as unknown as StaticModel;
    expect(shownModelIssues(model, () => true)).toEqual([]);
    expect(shownModelIssues(model, () => false)).toEqual([issue]);
  });

  it('`seq.resolve.system.unnamed` is said once the source opened', async () => {
    frames = installFrames();
    const v = createSeqViewer({ ...DefaultSeqViewerSpec(), ui: undefined });
    viewers.push(v);
    await v.load(spec('ring'));
    v.resize(600);
    await v.settle();
    const said = modelIssues(v);
    expect(codes(said)).toContain('seq.resolve.system.unnamed');
    // the wrapped row's length is never known: said too, its source settled
    expect(codes(said)).toContain('seq.resolve.length.unknown');
  });
});
