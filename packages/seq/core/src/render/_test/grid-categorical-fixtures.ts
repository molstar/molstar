/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// A grid on categorical axes (EC-35): an alignment's rows as both axes, a kinship-like grid over them.

import { PCG } from '@molstar/core/data/util/hash-functions';
import type { Node, Value } from '@molstar/core/repr/tree';
import { createSeqSpec } from '@molstar/seqx-core/lang';

/** Rows `S0…` with a random `score` each, a kinship-like grid over rows × rows (`i·n + j`), the plot in an order. */
export function kinship(
  n: number,
  seed: number,
  order: Value = { field: 'score' },
  orderY: Value = order,
  o: { readonly tree?: string } = {},
): { spec: Node; scores: number[] } {
  const r = new PCG(seed);
  const scores = Array.from({ length: n }, () => r.int() % 1000);
  const b = createSeqSpec();
  const labels = Array.from({ length: n }, (_, i) => `S${i}`);
  const aln = b.data.alignment('aln', {
    rows: { label: labels, letters: labels.map(() => 'ACGT'), score: scores },
  } as never);
  if (o.tree) b.data.tree('tree', { of: aln, newick: o.tree } as never);
  const rows = `${aln}/(rows)`;
  const values = new Float32Array(n * n);
  for (let i = 0; i < n; i++) for (let j = 0; j < n; j++) values[i * n + j] = i * n + j;
  const g = b.data.grid('kin', { x: rows, y: rows, values } as never);
  b.views.view('v', (v) => {
    v.section('plot', { layout: 'grid', of: rows, of_y: rows, member_order: order, member_order_y: orderY } as never)
      .track('kinship')
      .layer('g', { type: 'grid', data: g } as never);
  });
  return { spec: b.build(), scores };
}
