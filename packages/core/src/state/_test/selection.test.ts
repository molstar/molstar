/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { afterEach, beforeEach, describe, expect, it, vi } from 'vitest';
import { State, StateObject, StateSelection, StateTransform, StateTransformer } from '../index.js';
import { Task } from '@molstar/core/task';
import { ParamDefinition as PD } from '@molstar/core/util/param-definition';

interface TypeInfo {
  name: string;
  typeClass: 'Root' | 'Data';
}
const Create = StateObject.factory<TypeInfo>();
class Root extends Create({ name: 'Root', typeClass: 'Root' }) {}
class Leaf extends Create<{ value: number }>({ name: 'Leaf', typeClass: 'Data' }) {}

/** value 0: a Null object */
const Node = StateTransformer.create<Root | Leaf, Leaf, { value: number }>('selection-test', {
  name: 'node',
  from: [Root, Leaf],
  to: [Leaf],
  display: { name: 'Node' },
  params: () => ({ value: PD.Numeric(1) }),
  apply({ params }) {
    if (params.value === 0) return StateObject.Null;
    return new Leaf({ value: params.value });
  },
});

function newState() {
  return State.create(new Root({}), { runTask: <T>(t: Task<T>) => t.run() });
}

async function commit(state: State, edit: (b: ReturnType<State['build']>) => void) {
  const b = state.build();
  edit(b);
  await state.runTask(state.updateTree(b, { doNotLogTiming: true }));
}

describe('StateSelection', () => {
  beforeEach(() => {
    vi.spyOn(console, 'error').mockImplementation(() => {});
  });
  afterEach(() => {
    vi.restoreAllMocks();
  });

  it('subtree and children skip nodes that have no cell (perf-audit-2 03/B9)', async () => {
    const state = newState();
    await commit(state, (b) => {
      b.toRoot().apply(Node, { value: 0 }, { ref: 'null' });
    });
    // two edits: the nodes below the Null object get no cell
    await commit(state, (b) => {
      b.to('null').apply(Node, { value: 1 }, { ref: 'a', tags: 'x' });
      b.to('null').apply(Node, { value: 2 }, { ref: 'b', tags: 'x' });
    });
    expect(state.cells.has('a')).toBe(false);

    const subtree = state.select(StateSelection.Generators.byRef('null').subtree().withTag('x'));
    expect(subtree).toEqual([]);
    const all = state.select(StateSelection.Generators.byRef('null').subtree());
    expect(all.map((c) => c.transform.ref)).toEqual(['null']);
    const children = state.select(StateSelection.Generators.byRef('null').children().withTag('x'));
    expect(children).toEqual([]);
  });

  it('findTagInSubtree returns the first match in pre-order and stops there (perf-audit-2 03/B10)', async () => {
    const state = newState();
    await commit(state, (b) => {
      const first = b.toRoot().apply(Node, { value: 1 }, { ref: 'first', tags: 't' });
      first.apply(Node, { value: 1 }, { ref: 'first-child', tags: 't' });
      b.toRoot().apply(Node, { value: 1 }, { ref: 'second', tags: 't' });
      for (let i = 0; i < 50; i++) b.toRoot().apply(Node, { value: 1 }, { ref: `other-${i}` });
    });

    const get = vi.spyOn(state.tree.transforms, 'get');
    expect(StateSelection.findTagInSubtree(state.tree, StateTransform.RootRef, 't')).toBe('first');
    // the root and 'first'
    expect(get.mock.calls.length).toBe(2);
    expect(StateSelection.findTagInSubtree(state.tree, 'second', 't')).toBe('second');
    expect(StateSelection.findTagInSubtree(state.tree, StateTransform.RootRef, 'none')).toBeUndefined();
  });
});
