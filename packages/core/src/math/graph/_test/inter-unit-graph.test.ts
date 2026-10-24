/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { describe, expect, it } from 'vitest';
import { InterUnitGraph } from '../inter-unit-graph.js';

describe('InterUnitGraph', () => {
  it('getEdgeIndices of a vertex without edges allocates nothing (perf-audit-2 02/P9)', () => {
    const builder = new InterUnitGraph.Builder<number, number, { order: number }>();
    builder.startUnitPair(0, 1);
    builder.add(3, 5, { order: 1 });
    builder.finishUnitPair();
    const graph = new InterUnitGraph(builder.getMap());
    expect(graph.getEdgeIndices(3, 0)).toEqual([0]);
    const none = graph.getEdgeIndices(4, 0);
    expect(none).toEqual([]);
    expect(graph.getEdgeIndices(7, 1)).toBe(none);
  });
});
