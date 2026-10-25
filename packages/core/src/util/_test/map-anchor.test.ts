/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { setFlagsFromString } from 'node:v8';
import { runInNewContext } from 'node:vm';
import { describe, expect, it } from 'vitest';
import { anchorMap } from '../map-anchor.js';

describe('anchorMap', () => {
  it('keeps what it is given through a full GC; an object not given is collected (perf-audit-3 08/P12)', async () => {
    const [kept, dropped] = (() => {
      const value = { a: 1, b: 0.5 };
      expect(anchorMap(value)).toBe(value);
      return [new WeakRef(value), new WeakRef({ a: 1, b: 0.5 })];
    })();
    setFlagsFromString('--expose-gc');
    const gc: () => void = runInNewContext('gc');
    // WeakRef targets are kept until the job that made them ends
    await new Promise((r) => setTimeout(r, 0));
    gc();
    expect(kept.deref()).toEqual({ a: 1, b: 0.5 });
    expect(dropped.deref()).toBeUndefined();
  });
});
