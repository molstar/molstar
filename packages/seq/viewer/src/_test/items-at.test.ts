/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { describe, expect, it } from 'vitest';
import { blockUnder } from '../card-providers.js';
import { Describer, type DescribeSource } from '../describe.js';

// The items a card lists under a position (`Describer.itemsAt`): one item seen through two layers of the same data is
// one row; two distinct items that share a label and an extent are two. A system in two blocks: the block the position
// is on (`blockUnder`).

const item = (key: string) => ({ key, label: 'f', start: 10, end: 20, x0: 10, x1: 20 });

/** A layer of a data node (`address`) with its items. */
const layer = (key: string, address: string, items: ReturnType<typeof item>[]) =>
  ({ key, layer: { data: { address, system: 's' } }, type: { items: () => items }, lay: {} }) as never;

function section(layers: unknown[]) {
  return {
    axis: { x: (p: number) => p },
    section: { axis: { system: 's' } },
    tracks: [{ layers }],
    placement: { strips: { n: 0 } },
    lays: [],
  } as never;
}

const describer = () =>
  (
    new Describer({
      layout: () => ({}) as never,
      view: { model: { value: { systems: new Map() } } },
    } as unknown as DescribeSource) as unknown as { _words: { itemsAt: (s: never, p: number) => unknown[] } }
  )._words;

describe('itemsAt', () => {
  it('the same item through two layers of one data node is one row', () => {
    const a = layer('/a', '/data/bands', [item('q1')]);
    const b = layer('/b', '/data/bands', [item('q1')]);
    expect(describer().itemsAt(section([a, b]), 15)).toHaveLength(1);
  });

  it('two distinct items with one label and extent stay two: other data, or other keys', () => {
    const a = layer('/a', '/data/one', [item('x')]);
    const b = layer('/b', '/data/two', [item('x')]);
    expect(describer().itemsAt(section([a, b]), 15)).toHaveLength(2);
    const c = layer('/c', '/data/one', [item('x'), item('y')]);
    expect(describer().itemsAt(section([c]), 15)).toHaveLength(2);
  });
});

describe('blockUnder', () => {
  const block = (system: string, ...ranges: number[]) => ({ axis: { system, ranges } });

  it('one system in two blocks: the block whose range holds the position', () => {
    const blocks = [block('chr1', 0, 100), block('chr2', 0, 50), block('chr1', 200, 300)];
    expect(blockUnder(blocks, 'chr1', 50)).toBe(0);
    expect(blockUnder(blocks, 'chr1', 250)).toBe(2);
    // on none of its blocks: the first
    expect(blockUnder(blocks, 'chr1', 150)).toBe(0);
    expect(blockUnder(blocks, 'chr3', 5)).toBe(-1);
  });
});
