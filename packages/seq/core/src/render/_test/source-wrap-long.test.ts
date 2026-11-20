/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The row reader's wrap past one tile and past one lap (files.md, The file's own syntax): the origin's piece of a row
// ending past the length is read though the view is at the first positions of a system longer than a tile; a row longer
// than one lap is one lap; a start past the length wraps on a cyclic system and is left out on a linear one, said.

import type { Node } from '@molstar/core/repr/tree';
import { RecordsRuntime } from '@molstar/seqx-core/data/records';
import { createSeqSpec } from '@molstar/seqx-core/lang';
import { describe, expect, it } from 'vitest';
import { categoryColumn } from '@molstar/core/repr/tree';
// the render entry first: the items and ranges modules import each other
import '../index.js';
import type { RangesLay } from '../index.js';
import { collectItems, newItemRows } from '../items.js';
import { settle } from './helpers.js';

const Long = 3000;
const Size = { width: 800, height: 120 };
const row = (id: string, s: number, e: number) => `ring\tt\tCDS\t${s}\t${e}\t.\t+\t.\tID=${id}`;

function spec(gff: string, cyclic: boolean, data = (features: string) => features, length = Long): Node {
  const b = createSeqSpec();
  const ring = b.data.sequence('ring', { alphabet: 'dna', letters: 'ACGT'.repeat(length / 4), cyclic });
  const file = b.data.source('file', { format: 'gff3', value: gff });
  const features = b.data.ranges('features', { target: ring, source: file });
  b.views.view('main', (v) =>
    v
      .section('s', { of: ring })
      .track('t', { label: 'F' })
      .layer('l', { type: 'ranges', data: data(features) }),
  );
  return b.build();
}

async function laid(node: Node, viewport: { start: number; end: number }) {
  const s = await settle(node, [RecordsRuntime]);
  const f = await s.frame({ size: Size, settled: true, viewports: { '/views/main/s': viewport } });
  return { f, lay: f.layout.sections[0].tracks[0].layers[0].lay as RangesLay };
}

const pieces = (lay: RangesLay) =>
  Array.from({ length: lay.items }, (_, i) => `${lay.keys[i]}:${lay.starts[i]}-${lay.ends[i]}`).sort();

describe('a source’s rows past the end of a long system', () => {
  it('the origin’s piece of a row ending past the length is read from a view of the first positions', async () => {
    const { lay } = await laid(spec(row('across', 2901, 3100), true), { start: 0, end: 50 });
    expect(pieces(lay)).toContain('across:0-100');
  });

  it('a row longer than one lap is one lap, counted', async () => {
    const { f, lay } = await laid(spec(row('long', 11, 7000), true), { start: 0, end: Long });
    expect(pieces(lay)).toEqual(['long:0-10', 'long:10-3000']);
    const issue = f.layout.issues.find((i) => i.code === 'seq.render.item.range')!;
    expect(issue.message).toContain('1 longer than one lap');
  });

  it('a start at or past the length wraps on a cyclic system and is left out on a linear one, said (the reader’s rule)', () => {
    // no tile of a system holds such a row (a want is clipped to the length): the rule is the reader's for a source
    // that serves one, so it is read here from a table as such a source would give it
    const entry = {
      status: 'loaded',
      region: { source: 's', system: 'ring', level: 1, tile: 0, start: 0, end: Long },
      tables: {
        ranges: {
          form: 'list',
          start: Int32Array.of(3100, 2950),
          end: Int32Array.of(3200, 3100),
          id: categoryColumn(['past', 'over']),
        },
      },
    };
    const read = (cyclic: boolean) => {
      const rows = newItemRows();
      collectItems(rows, [entry as never], ['category'], undefined, undefined, { length: Long, cyclic });
      const list = Array.from({ length: rows.n }, (_, i) => `${rows.keys[i]}:${rows.starts[i]}-${rows.ends[i]}`).sort();
      return { list, rows };
    };
    expect(read(true).list).toEqual(['over:0-100', 'over:2950-3000', 'past:100-200']);
    const line = read(false);
    expect(line.list).toEqual(['over:2950-3000']);
    expect([line.rows.clipped, line.rows.leftOut]).toEqual([1, 1]);
  });

  it('data forms that are neither a base nor a trailing filter are said, every row drawn', async () => {
    const gff = [row('a', 1, 10), row('b', 20, 30)].join('\n');
    const node = spec(gff, false, (f) => `${f}/(except :of (item :id a))/(item :id b)`);
    const { f, lay } = await laid(node, { start: 0, end: Long });
    expect(f.layout.issues.map((i) => i.code)).toContain('seq.render.filter');
    expect(pieces(lay)).toEqual(['a:0-10', 'b:19-30']);
  });
});
