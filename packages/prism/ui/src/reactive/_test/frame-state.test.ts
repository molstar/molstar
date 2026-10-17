/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { describe, expect, it } from 'vitest';
import { bind, createRoot, flush, frameState, Scope, runInScope, watch } from '../../dom/scope.js';
import { flushFrame } from '../frame-state.js';
import { State } from '../state.js';

// The frame lane (rendering.md §2.3): offers keep the latest value, published once a frame; a watch following a frame
// state renders once in the frame, after the microtask lane; a paused scope skips both; a disposed offer is dropped.

describe('frame states', () => {
  it('ten offers in a frame: one publish, one render', () => {
    let renders = 0;
    let told = 0;
    const dispose = createRoot((dispose) => {
      const pointer = frameState(0);
      pointer.listen(() => told++);
      watch([pointer], () => renders++);
      for (let i = 1; i <= 10; i++) pointer.offer(i);
      expect(pointer.value).toBe(0);
      expect(pointer.latest).toBe(10);
      flush();
      expect(renders).toBe(1);
      flushFrame();
      expect(pointer.value).toBe(10);
      return dispose;
    });
    expect([told, renders]).toEqual([1, 2]);
    dispose();
  });

  it('catches a listener that is not coalesced (the check can fail): a state set ten times renders ten', () => {
    let writes = 0;
    const dispose = createRoot((dispose) => {
      const plain = new State(0);
      bind(plain, () => writes++);
      for (let i = 1; i <= 10; i++) plain.set(i);
      return dispose;
    });
    expect(writes).toBe(11);
    dispose();
  });

  it('a watch of a frame state and a model state renders once in the frame, after both', () => {
    let renders = 0;
    const model = new State('a');
    const dispose = createRoot((dispose) => {
      const pointer = frameState(0);
      watch([model, pointer], () => renders++);
      renders = 0;
      model.set('b');
      pointer.offer(1);
      flush();
      expect(renders).toBe(0);
      flushFrame();
      expect(renders).toBe(1);
      return dispose;
    });
    dispose();
  });

  it('a paused scope renders nothing in either lane; one render on resume', () => {
    let renders = 0;
    const root = new Scope();
    const pointer = runInScope(root, () => frameState(0));
    runInScope(root, () => watch([pointer], () => renders++));
    renders = 0;
    root.setActive(false);
    pointer.offer(5);
    flushFrame();
    expect(renders).toBe(0);
    root.setActive(true);
    expect(renders).toBe(1);
    root.dispose();
  });

  it('its scope gone, a pending offer is dropped', () => {
    const root = new Scope();
    const pointer = runInScope(root, () => frameState(0));
    pointer.offer(3);
    root.dispose();
    flushFrame();
    expect(pointer.value).toBe(0);
  });
});
