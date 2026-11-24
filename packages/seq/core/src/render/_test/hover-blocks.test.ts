/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// A hover over blocks outlines its position only on the strips of its own system (a CDK2 position lights both CDK2
// copies, never a cyclin whose axis holds the same number), canvas and SVG. (A ring of blocks draws no letters: its arcs carry the items, whose hover is keyed by element.)

import { describe, expect, it } from 'vitest';
import { type DrawLists, drawListsToSvg, recordDrawLists } from '../index.js';
import { settle } from './helpers.js';
import { blocksSpec } from './rings-fixtures.js';

const Size = { width: 900, height: 600 };
const Cdk2 = '/data/cdk2';

function hoverRects(lists: DrawLists): string[] {
  return recordDrawLists(lists, ['hover'])
    .split('\n')
    .filter((l) => l.startsWith('rect '));
}

const cases = [
  { name: 'each', options: { copies: 'each' as const } },
  { name: 'row', options: {} },
  { name: 'stack', options: { arrange: 'stack' as const, copies: 'each' as const } },
];

describe('the hover over blocks', () => {
  for (const c of cases) {
    for (const position of [10, 100]) {
      it(`${c.name}: CDK2 position ${position} lights the CDK2 copies, no cyclin`, async () => {
        const s = await settle(blocksSpec(c.options));
        const input = { size: Size };
        await s.frame(input);
        const hover = { layer: '/views/v/complex/letters/letters', section: 'complex', system: Cdk2, position };
        const f = s.renderer.frame({ ...input, emphasis: { hover } });
        // merged copies are one strip: one outline; each, two (a cyclin's would be a third)
        expect(hoverRects(f.lists).length).toBe(c.options.copies === 'each' ? 2 : 1);
        expect((drawListsToSvg(f.lists).match(/stroke="#2d271e/g) ?? []).length).toBe(
          c.options.copies === 'each' ? 2 : 1,
        );
      });
    }
  }
});
