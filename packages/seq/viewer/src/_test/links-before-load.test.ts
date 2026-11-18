/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { createSeqSpec } from '@molstar/seqx-core/lang';
import { afterEach, describe, expect, it } from 'vitest';
import { viewerItem } from '../config.js';
import { DefaultSeqViewerSpec } from '../default-spec.js';
import { createSeqViewer, type SeqViewer } from '../viewer.js';
import { type Frames, installFrames } from './frames.js';

// A link followed before the spec is loaded is kept, not refused (_spec/README.md "AddressModel"): it runs when the spec
// arrives, is kept again until the view is placed, and records its address as the focus when the kept move lands.

const Protein = '/views/main/protein';

function spec() {
  const b = createSeqSpec({ title: 'early' });
  const p = b.data.sequence('p', { alphabet: 'protein', letters: 'M'.repeat(40) });
  b.refs.ref('b', { to: `${p}/(range :from 10 :to 20)` });
  b.views.view('main', (v) => {
    v.section('protein', { of: p }).track('letters').layer('l', { type: 'letters' });
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

function make() {
  frames = installFrames();
  const viewer = createSeqViewer({
    ...DefaultSeqViewerSpec(),
    ui: undefined,
    config: [[viewerItem('viewport.maxPxPerPosition'), 1000]],
  });
  viewers.push(viewer);
  return viewer;
}

describe('a link followed before the spec is loaded', () => {
  it('waits, runs when the spec arrives, and records its focus when the view is placed', async () => {
    const v = make();
    const run = v.runLink('!seq-focus=/refs/b');
    expect(run?.outcomes.map((o) => [o.done, o.status])).toEqual([[false, 'waiting']]);
    await v.load(spec());
    // loaded, not placed: the move is kept as the view's intent, and says so
    expect(v.view.get(Protein)?.pending).toBe(true);
    expect(v.addresses.focusAddress()).not.toBe('/refs/b');
    v.resize(600);
    await v.settle();
    frames!.flush();
    expect(v.view.get(Protein)?.pending).toBe(false);
    expect([v.view.get(Protein)!.start, v.view.get(Protein)!.end]).toEqual([9, 20]);
    expect(v.addresses.focusAddress()).toBe('/refs/b');
  });

  it('a later move by hand ends the record: the focus reads as the shown region again', async () => {
    const v = make();
    v.runLink('!seq-focus=/refs/b');
    await v.load(spec());
    v.resize(600);
    await v.settle();
    v.view.pan(Protein, 5);
    expect(v.addresses.focusAddress()).not.toBe('/refs/b');
  });

  it('an address that names nothing in the spec loaded is refused', async () => {
    const v = make();
    expect(v.runLink('!seq-focus=/refs/nope')?.outcomes[0].status).toBe('waiting');
    await v.load(spec());
    v.resize(600);
    await v.settle();
    expect(v.addresses.resolve('/refs/nope').status).toBe('missing');
    expect(v.addresses.focus('/refs/nope')).toMatchObject({ done: false, message: 'Nothing to focus' });
  });
});
