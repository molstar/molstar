/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The tree's strokes lie inside the tree column (_spec/README.md "Render", The rows layer, The tree's look): the root's
// vertical, centred on the column's left edge, was clipped to half its width by the header's clip.

import { createSeqSpec } from '@molstar/seqx-core/lang';
import { describe, expect, it } from 'vitest';
import { recordDrawLists } from '../index.js';
import type { RowsLay } from '../rows.js';
import { settle } from './helpers.js';

const Labels = ['HBA/1-6', 'HBB/2-7', 'MYG/1-6', 'GLB/3-8', 'CYG/1-6'];
const Letters = ['MVLSPA', 'MVHLTP', 'MGLS--', 'MKAVLG', 'MEKVPG'];
const Newick = '((HBA:1,HBB:1)hb:1,(MYG:1,(GLB:1,CYG:1)gc:1)rest:1)root;';

function spec() {
  const b = createSeqSpec();
  const aln = b.data.alignment('aln', { alphabet: 'protein', rows: { label: Labels, letters: Letters } });
  const tree = b.data.tree('tree', { of: aln, newick: Newick });
  b.views.view('v', (v) =>
    v.section('s', { of: aln } as never, (s) =>
      s.track('rows', { label: '' }, (t) => t.layer('r', { type: 'rows', data: aln, tree, order: 'tree' } as never)),
    ),
  );
  return b;
}

describe('the tree column', () => {
  for (const dpr of [1, 2]) {
    it(`every stroke lies inside the column (dpr ${dpr})`, async () => {
      const s = await settle(spec().build(), []);
      const f = await s.frame({
        size: { width: 600, height: 300 },
        dpr,
        settled: true,
        viewports: {
          '/views/v/s': { start: 0, end: 6, bands: { '/views/v/s/rows': { bandTop: 0, rowScale: 20 } as never } },
        },
      });
      let lay: RowsLay | undefined;
      for (const sec of f.layout.sections)
        for (const l of sec.lays) for (const x of l ?? []) if (x.type.name === 'rows') lay = x.lay as RowsLay;
      expect(lay).toBeDefined();
      // the tree's rects: left of the labels
      const rects = recordDrawLists(f.lists, ['labels'])
        .split('\n')
        .filter((l) => l.startsWith('rect '))
        .map((l) =>
          l
            .slice(5)
            .split(' ')
            .map(Number),
        )
        .filter((r) => r[0] < lay!.labelX);
      expect(rects.length).toBeGreaterThan(3);
      const x0 = lay!.treeX0;
      const x1 = lay!.treeX1;
      for (const r of rects) {
        expect(r[0]).toBeGreaterThanOrEqual(x0 - 1e-9);
        expect(r[0] + r[2]).toBeLessThanOrEqual(x1 + 1e-9);
      }
    });
  }
});
