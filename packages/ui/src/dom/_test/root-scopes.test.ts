/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { afterEach, describe, expect, it } from 'vitest';
import { ViewScope } from '../scope.js';

// Root scopes register with the registry a test setup puts on `globalThis` (`vitest.setup.ts` disposes what a file
// leaves), and leave it when disposed; without one, nothing is kept.

const Key = Symbol.for('molstar.test.rootScopes');
const globals = globalThis as Record<symbol, unknown>;
const setupRegistry = globals[Key];
afterEach(() => {
  globals[Key] = setupRegistry;
});

function fakeRegistry() {
  const live = new Set<object>();
  globals[Key] = { add: (s: object) => live.add(s), delete: (s: object) => live.delete(s) };
  return live;
}

describe('root scopes', () => {
  it('register while live, not their children', () => {
    const live = fakeRegistry();
    const s = new ViewScope();
    const c = s.child();
    expect([...live]).toEqual([s]);
    c.dispose();
    expect([...live]).toEqual([s]);
    s.dispose();
    expect(live.size).toBe(0);
    // disposed again: nothing more
    s.dispose();
    expect(live.size).toBe(0);
  });

  it('keep nothing without a registry (apps)', () => {
    globals[Key] = undefined;
    const s = new ViewScope();
    expect(() => s.dispose()).not.toThrow();
  });
});
