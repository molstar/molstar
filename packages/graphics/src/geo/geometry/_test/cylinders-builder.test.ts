/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { describe, expect, it, vi } from 'vitest';

// counts the generic ChunkedArray adders the builder calls (it takes them once, at module load)
const calls = vi.hoisted(() => ({ n: 0 }));
vi.mock('@molstar/core/data/util', async (importOriginal) => {
  const m = await importOriginal<typeof import('@molstar/core/data/util')>();
  const { add, add3 } = m.ChunkedArray;
  const counted = <F extends (...args: any[]) => any>(f: F) =>
    ((...args: Parameters<F>) => {
      calls.n++;
      return f(...args);
    }) as F;
  return { ...m, ChunkedArray: { ...m.ChunkedArray, add: counted(add), add3: counted(add3) } };
});

const { CylindersBuilder } = await import('../cylinders/cylinders-builder.js');

/** cylinder i: start (i, 0, 0), end (i, 1, 0), scale i / 10, caps by i % 4, colorMode i % 3, group 100 + i */
function build(count: number, builder: ReturnType<typeof CylindersBuilder.create>) {
  for (let i = 0; i < count; i++) builder.add(i, 0, 0, i, 1, 0, i / 10, (i & 1) > 0, (i & 2) > 0, i % 3, 100 + i);
  return builder.getCylinders();
}

function expectVertices(cylinders: ReturnType<typeof build>, count: number) {
  expect(cylinders.cylinderCount).toBe(count);
  const start = cylinders.startBuffer.ref.value;
  const end = cylinders.endBuffer.ref.value;
  const group = cylinders.groupBuffer.ref.value;
  const scale = cylinders.scaleBuffer.ref.value;
  const cap = cylinders.capBuffer.ref.value;
  const colorMode = cylinders.colorModeBuffer.ref.value;
  for (let i = 0; i < count; i++) {
    for (let v = 6 * i; v < 6 * i + 6; v++) {
      expect([start[3 * v], start[3 * v + 1], start[3 * v + 2]]).toEqual([i, 0, 0]);
      expect([end[3 * v], end[3 * v + 1], end[3 * v + 2]]).toEqual([i, 1, 0]);
      expect(group[v]).toBe(100 + i);
      expect(scale[v]).toBe(Math.fround(i / 10));
      expect(cap[v]).toBe((i & 1 ? 1 : 0) + (i & 2 ? 2 : 0));
      expect(colorMode[v]).toBe(i % 3);
    }
  }
}

describe('CylindersBuilder.add', () => {
  it('writes the 6 vertices of a cylinder into the current chunks, not through the adders (perf-audit-2 08/P16)', () => {
    calls.n = 0;
    const cylinders = build(100, CylindersBuilder.create(100, 50));
    expect(calls.n).toBe(0);
    expectVertices(cylinders, 100);
  });

  it('crosses chunk ends, and refills reused buffers', () => {
    // chunks of 2 cylinders after a first of 1: the adders allocate each next one
    const cylinders = build(7, CylindersBuilder.create(1, 2));
    expectVertices(cylinders, 7);
    expectVertices(build(5, CylindersBuilder.create(5, 2, cylinders)), 5);
  });
});
