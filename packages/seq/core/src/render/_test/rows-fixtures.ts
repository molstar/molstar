/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { PCG } from '@molstar/core/data/util/hash-functions';
import type { Node } from '@molstar/core/repr/tree';
import { cellsSource, type FormatRuntime } from '@molstar/seqx-core/data';
import { createSeqSpec } from '@molstar/seqx-core/lang';

// A synthetic alignment of any number of rows (300 columns), its rows made by index, and a spec of its rows layer: the
// rows' time and allocation cases share it.

const Residues = 'ACDEFGHIKLMNPQRSTVWY-';
export const Columns = 300;
export const Track = '/views/v/s/rows';
export const Section = '/views/v/s';

/** A synthetic alignment source (`synthetic`): rows made by index, each a biased draw a column. */
export const Synthetic: FormatRuntime = {
  name: 'synthetic',
  formats: ['synthetic'],
  async open(source) {
    const n = (source.value as { rows: number }).rows;
    return cellsSource(n, Columns, (start, end) => {
      const codes = new Uint8Array((end - start) * Columns);
      for (let i = start; i < end; i++) {
        const rng = new PCG((i * 2654435761) >>> 0);
        for (let j = 0; j < Columns; j++)
          codes[(i - start) * Columns + j] = Residues.charCodeAt(
            rng.float() < 0.6 ? (j * 7) % Residues.length : rng.int() % Residues.length,
          );
      }
      const label = Array.from({ length: end - start }, (_, k) => `R${start + k}`);
      return { codes, label, first: new Int32Array(end - start).fill(1) };
    });
  },
};

/** A balanced tree over the rows, its leaves `R<i>`. */
function balanced(n: number): string {
  const out: string[] = [];
  const grow = (a: number, b: number) => {
    if (b - a === 1) return void out.push(`R${a}:1`);
    const m = (a + b) >> 1;
    out.push('(');
    grow(a, m);
    out.push(',');
    grow(m, b);
    out.push('):1');
  };
  grow(0, n);
  return `${out.join('')};`;
}

export function rowsSpec(rows: number, tree: boolean): Node {
  const b = createSeqSpec();
  const src = b.data.source('src', { format: 'synthetic', value: { rows } });
  const aln = b.data.alignment('aln', { alphabet: 'protein', source: src });
  const t = tree ? b.data.tree('tree', { of: aln, newick: balanced(rows) }) : undefined;
  b.views.view('v', (v) =>
    v
      .section('s', { of: aln })
      .track('rows', { label: '', height: 600 })
      .layer('r', { type: 'rows', data: aln, ...(t ? { tree: t, order: 'tree' } : {}) } as never),
  );
  return b.build();
}
