/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// `useDerived` computes once, then again only when a state or one of its `deps` changed: not on every render, and not
// on each of React's snapshot reads of one render.

import { act } from 'react';
import { describe, expect, it } from 'vitest';
import { useDerived } from '../../src/react/index.js';
import { state } from '../../src/state.js';
import { harness } from './fixtures.js';

const { parent } = harness();

describe('useDerived', () => {
  it('computes once across renders for unrelated reasons, and once per change of a state or a dep', async () => {
    const [a, b] = [state(1), state(10)];
    let computes = 0;
    let renders = 0;
    function Part({ factor, unrelated }: { factor: number; unrelated: number }) {
      const sum = useDerived(
        [a, b],
        () => {
          computes++;
          return a.value + b.value * factor;
        },
        { deps: [factor] },
      );
      renders++;
      return <span data-unrelated={unrelated}>{sum}</span>;
    }
    let factor = 1;
    const t = parent((i) => <Part factor={factor} unrelated={i} />);
    await t.render();
    expect([computes, t.container.textContent]).toEqual([1, '11']);
    // Three renders for a prop the compute does not read.
    await t.render();
    await t.render();
    await t.render();
    expect(renders).toBeGreaterThanOrEqual(4);
    expect(computes).toBe(1);
    // A state changed: computed again once.
    await act(async () => a.set(2));
    expect([computes, t.container.textContent]).toEqual([2, '12']);
    // A dep changed.
    factor = 2;
    await t.render();
    expect([computes, t.container.textContent]).toEqual([3, '22']);
    await t.render();
    expect(computes).toBe(3);
  });

  it('sees a state that changed between a render and the subscription', async () => {
    const a = state(1);
    let computes = 0;
    function Part() {
      const v = useDerived([a], () => {
        computes++;
        return a.value;
      });
      // Changes after this render's snapshot, before its effects subscribe.
      if (computes === 1 && a.value === 1) a.set(2);
      return <span>{v}</span>;
    }
    const t = parent(() => <Part />);
    await t.render();
    expect(t.container.textContent).toBe('2');
  });
});
