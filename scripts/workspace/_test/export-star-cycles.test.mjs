/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { describe, expect, it } from 'vitest';
import { components, starsInCycles } from '../export-star-cycles.mjs';

const graph = (edges) => {
  const g = new Map();
  for (const [from, to] of edges) {
    if (!g.has(from)) g.set(from, new Set());
    if (!g.has(to)) g.set(to, new Set());
    g.get(from).add(to);
  }
  return g;
};

describe('export * inside import cycles', () => {
  it('finds the components of a graph', () => {
    const c = components(
      graph([
        ['a', 'b'],
        ['b', 'c'],
        ['c', 'a'],
        ['c', 'd'],
      ]),
    );
    expect(c.get('a')).toBe(c.get('b'));
    expect(c.get('b')).toBe(c.get('c'));
    expect(c.get('d')).not.toBe(c.get('a'));
  });

  it('reports a star export on a cycle, not one off it, nor one whose cycle runs through import()', () => {
    // index -> mmcif (star) -> format -> index: a cycle; index -> rings (star): none; lazy -> index only by import()
    const static_ = graph([
      ['index', 'mmcif'],
      ['mmcif', 'format'],
      ['format', 'index'],
      ['index', 'rings'],
      ['index', 'lazy'],
    ]);
    const stars = [
      { from: 'index', to: 'mmcif' },
      { from: 'index', to: 'rings' },
      { from: 'index', to: 'lazy' },
    ];
    expect(starsInCycles({ static: static_, stars })).toEqual([{ from: 'index', to: 'mmcif' }]);
  });
});
