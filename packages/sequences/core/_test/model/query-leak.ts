/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The query heap guards' positive control (the guard template, lib/seqstar/AGENTS.md, Tests): `(leak Q)`, an
// extension's runtime around a query that makes one object a run of its selection, through `leak` (garbage at once,
// never a kept list).

import { defineSymbol } from '@molstar/seq-core/query';
import type { QueryExtension, QueryFn, SeqSelection } from '../../src/model/query/index.js';
import { leak } from './heap.js';

export const leakSymbol = defineSymbol({
  ...{ name: 'test.leak', alias: 'leak', kind: 'modifier', result: 'query', signature: '(leak Q)' },
  ...{ args: [{ name: 'Q', type: 'query', positional: true, doc: 'the query' }], doc: 'one object per run' },
  ...{ examples: [], hints: [] },
});

/** The `test` extension holding `(leak Q)`, with the guard's own letter classes. */
export const leakExtension = (classes: Record<string, Record<string, string>> = {}): QueryExtension => ({
  namespace: 'test',
  symbols: [leakSymbol],
  classes,
  runtimes: [
    {
      symbol: leakSymbol,
      compile: (c) => {
        const q = c.query('0') as QueryFn<SeqSelection>;
        return (ctx) => {
          const s = q(ctx);
          for (let k = 0; k < s.runCount; k++) leak(k);
          return s;
        };
      },
    },
  ],
});
