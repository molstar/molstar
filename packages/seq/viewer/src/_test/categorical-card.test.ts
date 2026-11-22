// @vitest-environment jsdom
/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// A categorical grid's card (EC-35; _spec/README.md "Behaviors", hoverCard): a cell named by its two members' labels
// and its value, whatever order the axes show them in.

import { createSeqSpec } from '@molstar/seqx-core/lang';
import { DomHost } from '@molstar/ui/dom';
import { afterEach, beforeEach, describe, expect, it, vi } from 'vitest';
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

const Rows = '/data/aln/(rows)';

function kinshipSpec() {
  const b = createSeqSpec({ title: 'kinship' });
  const labels = ['NA12878', 'NA12891', 'NA12892'];
  b.data.alignment('aln', { rows: { label: labels, letters: labels.map(() => 'ACGT') } } as never);
  const values = Float32Array.from([1, 0.42, 0.5, 0.42, 1, 0.01, 0.5, 0.01, 1]);
  const g = b.data.grid('kin', { x: Rows, y: Rows, values } as never);
  b.views.view('v', (v) => {
    v.section('plot', {
      layout: 'grid',
      of: Rows,
      of_y: Rows,
      member_order: { field: 'label', direction: 'descending' },
    } as never)
      .track('kinship')
      .layer('g', { type: 'grid', data: g } as never);
  });
  return b.build();
}

describe('a categorical grid’s card (EC-35)', () => {
  it('names a cell by its two members and its value', async () => {
    const v = createSeqViewer(DefaultSeqViewerSpec());
    viewers.push(v);
    await v.load(kinshipSpec());
    v.resize(900);
    await v.settle();
    frames.flush();
    v.draw();
    const s = v.frame!.layout.sections[0];
    const l = s.lays[0][0];
    const owner = {};
    v.hover.set(
      {
        layer: l.key,
        section: s.section.address,
        system: Rows,
        position: 0,
        ySystem: Rows,
        yPosition: 1,
        target: 'cell',
      } as never,
      owner,
    );
    const words = new Describer({
      spec: v.spec,
      store: v.store,
      addresses: v.addresses,
      view: v.view,
      layout: () => v.frame?.layout,
    });
    const card = words.describe(v.hover.current)!;
    expect(card.title).toBe('NA12878 × NA12891');
    expect(card.rows).toContainEqual({ label: 'Value', value: '0.42' });
    v.hover.clear(owner);
  });
});
