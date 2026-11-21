/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { createSeqSpec, theme } from '@molstar/seqx-core/lang';
import { afterEach, describe, expect, it } from 'vitest';
import { viewerItem } from '../config.js';
import { DefaultSeqViewerSpec } from '../default-spec.js';
import { legendRows } from '../legends.js';
import { createSeqViewer, type SeqViewer } from '../viewer.js';
import { type Frames, installFrames } from './frames.js';

// A ramp chosen from a layer's Color menu takes the layer's field domain (the data's extent when none is declared),
// not [0, 1]: pLDDT on 0 to 100 (the author's look, 2026-11-21).

let frames: Frames | undefined;
const viewers: SeqViewer[] = [];
afterEach(() => {
  for (const v of viewers.splice(0)) v.retire();
  frames?.restore();
  frames = undefined;
});

const Layer = '/views/main/p/conf/v';

async function make() {
  frames = installFrames();
  const b = createSeqSpec();
  const p = b.data.sequence('p', { alphabet: 'protein', letters: 'MEEPQSDPSVEPPLSQETFSDLWKLL' });
  const conf = b.data.values('conf', { target: p, values: Float32Array.from({ length: 26 }, (_, i) => 20 + i * 3) });
  b.views.view('main', (v) => {
    const s = v.section('p', { of: p });
    s.track('conf', { label: 'pLDDT' }).layer('v', {
      type: 'values',
      data: conf,
      color: theme.field('value', {
        palette: theme.binned({ edges: [0, 50, 70, 90, 100], colors: ['#ff7d45', '#ffdb13', '#65cbf3', '#0053d6'] }),
      }),
    });
  });
  const v = createSeqViewer({
    ...DefaultSeqViewerSpec(),
    ui: undefined,
    config: [[viewerItem('viewport.maxPxPerPosition'), 1000]],
  });
  viewers.push(v);
  await v.load(b.build());
  v.resize(600);
  await v.settle();
  frames.flush();
  return v;
}

describe('a ramp chosen from the Color menu', () => {
  it('takes the layer’s field domain (the data’s extent), not 0 to 1', async () => {
    const v = await make();
    await v.execute('layer.channel', { address: Layer, scheme: 'viridis' });
    await v.settle();
    frames!.flush();
    const ramp = v.legends.legends.value.find((l) => l.entries.kind === 'ramp')!;
    // the data's extent (20 to 95), held once the layout has folded it
    expect(ramp.domain).toEqual([20, 95]);
  });
});
