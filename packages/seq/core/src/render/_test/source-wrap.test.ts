/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// A source's row ending past its system's length (files.md, The file's own syntax): two rows of one item on a cyclic
// system, drawn as one item across the origin; clipped on a linear one, and said. A layer's leading item form is the
// track's own: its rejected rows are not counted in the "n filtered" chip, a trailing `(not …)` is.

import type { Node } from '@molstar/core/repr/tree';
import { RecordsRuntime } from '@molstar/seqx-core/data/records';
import { createSeqSpec } from '@molstar/seqx-core/lang';
import { describe, expect, it } from 'vitest';
import { hitView, type RangesLay } from '../index.js';
import { settle } from './helpers.js';

const Length = 100;
const Gff = [
  '##gff-version 3',
  'ring\tt\tCDS\t91\t110\t.\t+\t.\tID=across',
  'ring\tt\tCDS\t21\t30\t.\t+\t.\tID=inside',
  'ring\tt\tpromoter\t41\t50\t.\t+\t.\tID=prom',
  'ring\tt\tmisc\t51\t60\t.\t+\t.\tID=misc',
].join('\n');

function spec(cyclic: boolean, query = ''): Node {
  const b = createSeqSpec();
  const letters = 'ACGT'.repeat(Length / 4);
  const ring = b.data.sequence('ring', { alphabet: 'dna', letters, cyclic });
  const file = b.data.source('file', { format: 'gff3', value: Gff, fields: { category: 'type' } });
  const features = b.data.ranges('features', { target: ring, source: file });
  b.views.view('main', (v) =>
    v
      .section('s', { of: ring })
      .track('t', { label: 'Features' })
      .layer('l', { type: 'ranges', data: `${features}${query}` }),
  );
  return b.build();
}

const Size = { width: 800, height: 120 };

async function layout(cyclic: boolean, query = '', viewport = { start: 0, end: Length }) {
  const s = await settle(spec(cyclic, query), [RecordsRuntime]);
  const f = await s.frame({ size: Size, settled: true, viewports: { '/views/main/s': viewport } });
  const t = f.layout.sections[0].tracks[0];
  return { f, t, lay: t.layers[0].lay as RangesLay };
}

const pieces = (lay: RangesLay) =>
  Array.from({ length: lay.items }, (_, i) => `${lay.keys[i]}:${lay.starts[i]}-${lay.ends[i]}`).sort();

describe('a row past the end of its system', () => {
  it('on a cyclic system is two rows of one item, whether the origin is in view or not', async () => {
    // the origin in view: the tail's piece comes from the tile of the system's last position
    const whole = await layout(true);
    expect(pieces(whole.lay)).toEqual(['across:0-10', 'across:90-100', 'inside:20-30', 'misc:50-60', 'prom:40-50']);
    // a window on the first positions only: the item's second piece is read though its row starts at the end
    const start = await layout(true, '', { start: 0, end: 15 });
    expect(pieces(start.lay)).toContain('across:0-10');
    expect(whole.f.layout.issues.map((i) => i.code)).not.toContain('seq.render.item.range');
  });

  it('is one item: both pieces answer to one key, a hit on either names it', async () => {
    const { f, lay } = await layout(true);
    const keys = new Set(Array.from({ length: lay.items }, (_, i) => String(lay.keys[i])));
    expect(keys.has('across')).toBe(true);
    const box = f.layout.sections[0].tracks[0].box;
    const y = box.y + box.height / 2;
    const at = (position: number) => {
      const axis = f.layout.sections[0].axis!;
      return hitView(f.layout, axis.xOfUnit(axis.unitOf(0, position + 0.5)), y);
    };
    for (const position of [95, 5]) {
      const hit = at(position);
      expect(hit?.element?.key, `${position}`).toBe('across');
    }
  });

  it('on a linear system is clipped to its length, and said once', async () => {
    const { f, lay } = await layout(false);
    expect(pieces(lay)).toEqual(['across:90-100', 'inside:20-30', 'misc:50-60', 'prom:40-50']);
    const issues = f.layout.issues.filter((i) => i.code === 'seq.render.item.range');
    expect(issues.length).toBe(1);
    expect(issues[0].message).toContain('1 clipped');
  });
});

describe('a leading item form and a trailing filter', () => {
  it('a leading item form is the track: its rejected rows are not counted as filtered', async () => {
    const { t, lay } = await layout(true, '/(item :category CDS)');
    expect(lay.keys.slice(0, lay.items).sort()).toEqual(['across', 'across', 'inside']);
    expect(t.filter).toBeUndefined();
  });

  it('a trailing `(not …)` after it is the user’s: only what it rejects is counted', async () => {
    const { t, lay } = await layout(true, '/(item :category CDS)/(except :of (item :id inside))');
    expect(lay.keys.slice(0, lay.items).sort()).toEqual(['across', 'across']);
    expect(t.filter).toMatchObject({ count: 1, text: '1 filtered' });
    // `(not …)` alone filters as before: every row it rejects is counted
    const alone = await layout(true, '/(except :of (item :category CDS))');
    expect(alone.t.filter).toMatchObject({ count: 2 });
  });
});
