// @vitest-environment jsdom
/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The hover card over data read by region (_spec/README.md "Hover card", By region): a matrix cell as its item, its
// state and its value, by the core's words; the card and the live region alike.

import { createSeqSpec } from '@molstar/seqx-core/lang';
import type { Hit as RenderHit } from '@molstar/seqx-core/render';
import { DomHost } from '@molstar/ui/dom';
import { afterEach, beforeEach, describe, expect, it, vi } from 'vitest';
import { Announcer } from '../announcer.js';
import { RecordingContext } from '../canvas/_test/recording.js';
import { DefaultSeqViewerSpec } from '../default-spec.js';
import { Describer } from '../describe.js';
import { createSeqViewer, type SeqViewer } from '../viewer.js';
import { type Frames, installFrames } from './frames.js';

let frames: Frames;
const viewers: SeqViewer[] = [];

beforeEach(() => {
  frames = installFrames();
  vi.spyOn(DomHost, 'now').mockImplementation(() => 1000);
  vi.spyOn(HTMLCanvasElement.prototype, 'getContext').mockImplementation(function (this: HTMLCanvasElement) {
    return ((this as never as { __ctx?: RecordingContext }).__ctx ??= new RecordingContext()) as never;
  });
});

afterEach(() => {
  for (const v of viewers.splice(0)) v.retire();
  frames.restore();
  vi.restoreAllMocks();
});

/** Methylation over 8 positions in two states, one value missing. */
function spec() {
  const b = createSeqSpec({ title: 'methylation' });
  const s = b.data.sequence('s', { alphabet: 'dna', letters: 'ACGTACGT' });
  const values = Float32Array.from([0, 1, 0.5, 0.5, 0.2, Number.NaN, 1, 0, 0.1, 0.9, 0.3, 0.3, 0.7, 0.2, 0.4, 0.6]);
  const m = b.data.matrix('m', { target: s, values, states: ['tumour', 'normal'] });
  b.views.view('main', (v) => {
    v.section('s', { of: s }).track('m', { height: 40 }, (t) => {
      t.layer('m', { type: 'matrix', data: m });
    });
  });
  return b.build();
}

async function make() {
  const v = createSeqViewer(DefaultSeqViewerSpec());
  viewers.push(v);
  await v.load(spec());
  v.resize(900);
  await v.settle();
  frames.flush();
  v.draw();
  return v;
}

function cardOf(v: SeqViewer, hit: Partial<RenderHit> & { system: string }) {
  const owner = {};
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

describe('the hover card over data read by region', () => {
  it('names a matrix cell by its item, its state and its value; a missing one says no data', async () => {
    const v = await make();
    const s = v.frame!.layout.sections[0];
    let key = '';
    for (const ls of s.lays) for (const l of ls) if (l.type.name === 'matrix') key = l.key;
    expect(key).not.toBe('');
    const at = { layer: key, section: s.section.address, system: '/data/s', target: 'cell' as const };
    const statesSystem = '/data/m/(states)';
    const { card, said } = cardOf(v, { ...at, position: 2, ySystem: statesSystem, yPosition: 0 });
    expect(card?.rows).toEqual([
      { label: 'Item', value: 's:3' },
      { label: 'State', value: 'tumour' },
      { label: 'Value', value: '0.2' },
    ]);
    expect(said).toBe(card!.text);
    const missing = cardOf(v, { ...at, position: 2, ySystem: statesSystem, yPosition: 1 }).card;
    expect(missing?.text).toMatch(/normal, no data$/);
  });
});
