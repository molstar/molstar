/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The headers' allocation target (_spec/README.md "Render", Track headers): a pointer moving down the label column of
// 200 tracks, each with two slots, the "⋯" and a state, hit-tested and its emphasis drawn (the actions of the track it
// is over shown, an emphasis redraw), allocates nothing and lays nothing out. By hand, never in CI:
// `pnpm -C packages/seq/core test:perf`.

import { createSeqSpec } from '@molstar/seqx-core/lang';
import { quickCase } from '@molstar/test-harness/quick';
import { describe, expect, inject, it } from 'vitest';
import { type Emphasis, type HeaderPartsOf, hitView } from '../index.js';
import { P53 } from './fixtures.js';
import { settle } from './helpers.js';

const Tracks = 200;
const Cross = { commands: [0, 1, 0, 1], points: [4, 4, 20, 20, 4, 20, 20, 4] };
const Slots = [
  { name: 'seq:hide', label: 'Hide track', icon: Cross, enabled: true },
  { name: 'seq:color', label: 'Color', icon: Cross, enabled: false },
];
const States = [{ name: 'host:bound', text: 'bound' }];
const Host: HeaderPartsOf = (f) => (f.track ? { slots: Slots, states: States, menu: true } : undefined);

function spec() {
  const b = createSeqSpec();
  const p = b.data.sequence('p', { alphabet: 'protein', letters: P53 });
  b.views.view('main', (v) => {
    const s = v.section('p', { of: p });
    for (let i = 0; i < Tracks; i++) s.track(`t${i}`, { label: `Track ${i}` }).layer('l', { type: 'letters', data: p });
  });
  return b.build();
}

describe('the headers allocate nothing on a pointer move', () => {
  it('a pointer down the label column of 200 tracks: hit, the actions of its track drawn', async ({ task }) => {
    const s = await settle(spec(), [], { labelColumn: { min: 140 } });
    const size = { width: 1200, height: Tracks * 24 + 40 };
    const emphasis: { hover?: Emphasis['hover']; header?: string } = {};
    const input = { size, settled: true, headers: Host, emphasis };
    const f = await s.frame(input);
    const tracks = f.layout.sections[0].tracks;
    // the points (made before the walk): across each track's line, its label, its slots and its state
    const xs = [12, 40, ...tracks[0].header!.slots.map((x) => x.x + 4), tracks[0].header!.menu!.x + 4];
    const points: number[] = [];
    for (const t of tracks) for (const x of xs) points.push(x, t.well.y + t.well.height / 2);
    const addresses = tracks.map((t) => t.track.address);
    const n = points.length / 2;
    const step = (i: number) => {
      const k = i % n;
      emphasis.hover = hitView(f.layout, points[2 * k], points[2 * k + 1]);
      // the track the pointer is over (the viewer's headers say it from the geometry)
      emphasis.header = addresses[(k / xs.length) | 0];
      s.renderer.frame(input);
    };
    for (let i = 0; i < n; i++) {
      step(i);
      expect(s.renderer.frame(input).rebuilt).toBe(false);
    }
    const { r } = quickCase(task, 'headers pointer move', step, {});
    console.info(`headers pointer move: ${r.bytesPerStep.toFixed(0)} bytes a step (${r.note})`);
    expect(r.bytesPerStep, `headers pointer move: ${r.note}`).toBeLessThanOrEqual(inject('perfThreshold'));
  });
});
