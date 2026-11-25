// @vitest-environment jsdom
/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// A click on a matrix cell selects its item as a ranges click does (_spec/README.md "Behaviors"): the matrix's y is
// its states, not a grid's second axis, so the click is no grid cell.

import { createSeqSpec } from '@molstar/seqx-core/lang';
import { partOf } from '@molstar/seqx-core/query';
import type { MatrixLay } from '@molstar/seqx-core/render';
import { DomHost } from '@molstar/ui/dom';
import { afterEach, beforeEach, describe, expect, it, vi } from 'vitest';
import { RecordingContext } from '../canvas/_test/recording.js';
import { DefaultSeqViewerSpec } from '../default-spec.js';
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

async function make() {
  const b = createSeqSpec({ title: 'methylation' });
  const s = b.data.sequence('s', { alphabet: 'dna', letters: 'ACGTACGT' });
  const values = Float32Array.from([0, 1, 0.5, 0.5, 0.2, 0.1, 1, 0, 0.1, 0.9, 0.3, 0.3, 0.7, 0.2, 0.4, 0.6]);
  const m = b.data.matrix('m', { target: s, values, states: ['tumour', 'normal'] });
  b.views.view('main', (v) => {
    v.section('s', { of: s }).track('m', { height: 40 }, (t) => {
      t.layer('m', { type: 'matrix', data: m });
    });
  });
  const v = createSeqViewer(DefaultSeqViewerSpec());
  viewers.push(v);
  await v.load(b.build());
  v.resize(900);
  await v.settle();
  frames.flush();
  v.draw();
  return v;
}

const clickAt = (v: SeqViewer, x: number, y: number) => {
  v.input.move(x, y);
  v.input.down(x, y, {});
  v.input.up();
  v.input.click(x, y, {});
  frames.flush();
};

describe('a click on a matrix cell', () => {
  it('selects its item: its range on the data’s system, no grid cell', async () => {
    const v = await make();
    let lay: MatrixLay | undefined;
    for (const ls of v.frame!.layout.sections[0].lays)
      for (const l of ls) if (l.type.name === 'matrix') lay = l.lay as MatrixLay;
    expect(lay).toBeDefined();
    const a = lay!.axis!;
    const x = (a.x(2) + a.x(3)) / 2;
    clickAt(v, x, lay!.box.y + lay!.rowScale * 1.5 - lay!.bandTop);
    const sel = v.selection.value.value;
    expect(sel.location.rects).toEqual([]);
    expect([...(partOf(sel.location, '/data/s')?.ranges ?? [])]).toEqual([2, 3]);
  });
});
