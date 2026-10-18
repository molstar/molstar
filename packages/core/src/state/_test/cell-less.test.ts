/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { afterEach, beforeEach, describe, expect, it, vi } from 'vitest';
import { State, StateObject, StateTransformer } from '../index.js';
import { Task } from '@molstar/core/task';
import { ParamDefinition as PD } from '@molstar/core/util/param-definition';

// Nothing below a Null object is evaluated, and a node added below one in a multi-edit builder gets no cell (the diff
// does not descend there). Such nodes are reached anyway through dependencies: a Null node re-evaluated because a
// dependency changed, or a node that depends on one of them.

interface TypeInfo {
  name: string;
  typeClass: 'Root' | 'Data';
}
const Create = StateObject.factory<TypeInfo>();
class Root extends Create({ name: 'Root', typeClass: 'Root' }) {}
class Leaf extends Create<{ value: number }>({ name: 'Leaf', typeClass: 'Data' }) {}

/** value 0: a Null object; the value adds the parent's and the dependencies' values */
const Node = StateTransformer.create<Root | Leaf, Leaf, { value: number }>('cell-less-test', {
  name: 'node',
  from: [Root, Leaf],
  to: [Leaf],
  display: { name: 'Node' },
  params: () => ({ value: PD.Numeric(1) }),
  apply({ a, params, dependencies }) {
    if (params.value === 0) return StateObject.Null;
    let value = ((a.data as { value?: number }).value ?? 0) + params.value;
    for (const d of Object.values(dependencies ?? {})) value += (d as Leaf).data.value;
    return new Leaf({ value });
  },
});

function sumOf(dependencies: Record<string, StateObject> | undefined) {
  let value = 0;
  for (const d of Object.values(dependencies ?? {})) value += (d as Leaf).data.value;
  return value;
}

/** A Null object while its dependencies add up to 1 or less */
const Gate = StateTransformer.create<Root | Leaf, Leaf, {}>('cell-less-test', {
  name: 'gate',
  from: [Root, Leaf],
  to: [Leaf],
  display: { name: 'Gate' },
  params: () => ({}),
  apply({ dependencies }) {
    const value = sumOf(dependencies);
    return value <= 1 ? StateObject.Null : new Leaf({ value });
  },
});

/** Depends on `source` (a dependency derived from params, `getDependencies`); empty by default in the builder */
const Uses = StateTransformer.create<Root | Leaf, Leaf, { value: number; source: string }>('cell-less-test', {
  name: 'uses',
  from: [Root, Leaf],
  to: [Leaf],
  display: { name: 'Uses' },
  params: () => ({ value: PD.Numeric(1), source: PD.Text('') }),
  getDependencies: (params) => (params.source ? [params.source] : []),
  apply({ params, dependencies }) {
    return new Leaf({ value: params.value + sumOf(dependencies) });
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

const waiters = (state: State) =>
  (state as unknown as { dependencyWaiters: { byCell: Map<unknown, unknown> } }).dependencyWaiters;

const valueOf = (state: State, ref: string) => (state.cells.get(ref)?.obj?.data as { value?: number })?.value;

describe('State cells below Null objects', () => {
  // the state logs each transform error
  beforeEach(() => {
    vi.spyOn(console, 'error').mockImplementation(() => {});
  });
  afterEach(() => {
    vi.restoreAllMocks();
  });

  it('a Null node re-evaluated for a dependency has cells for its new children', async () => {
    const state = newState();
    await commit(state, (b) => {
      b.toRoot().apply(Node, { value: 1 }, { ref: 'dep' });
      b.toRoot().apply(Node, { value: 0 }, { ref: 'null', dependsOn: ['dep'] });
    });
    expect(state.cells.get('null')!.obj).toBe(StateObject.Null);

    // 'child' is below a Null object; 'null' is re-evaluated because 'dep' changed
    await commit(state, (b) => {
      b.to('null').apply(Node, { value: 1 }, { ref: 'child' });
      b.to('dep').update({ value: 2 });
    });
    expect(state.cells.get('null')!.obj).toBe(StateObject.Null);
    expect(state.cells.has('child')).toBe(true);
    expect(state.cells.get('child')!.obj).toBeUndefined();

    // once 'null' has an object, its child is evaluated
    await commit(state, (b) => {
      b.to('null').update({ value: 3 });
    });
    expect(valueOf(state, 'null')).toBe(5);
    expect(valueOf(state, 'child')).toBe(6);
  });

  it('a node that depends on a node without a cell is in error until that cell exists', async () => {
    const state = newState();
    await commit(state, (b) => {
      b.toRoot().apply(Node, { value: 0 }, { ref: 'null' });
    });
    // two edits: 'below' gets no cell
    await commit(state, (b) => {
      b.to('null').apply(Node, { value: 2 }, { ref: 'below' });
      b.toRoot().apply(Node, { value: 1 }, { ref: 'other' });
    });
    expect(state.cells.has('below')).toBe(false);

    await commit(state, (b) => {
      b.toRoot().apply(Node, { value: 1 }, { ref: 'user', dependsOn: ['below'] });
    });
    const user = state.cells.get('user')!;
    expect(user.status).toBe('error');
    expect(user.dependencies.dependsOn).toEqual([]);

    // 'below' is evaluated: 'user' links to it and is evaluated again
    await commit(state, (b) => {
      b.to('null').update({ value: 1 });
    });
    expect(valueOf(state, 'below')).toBe(3);
    expect(user.status).toBe('ok');
    expect(user.dependencies.dependsOn.map((c) => c.transform.ref)).toEqual(['below']);
    expect(state.cells.get('below')!.dependencies.dependentBy).toEqual([user]);
    expect(valueOf(state, 'user')).toBe(4);
  });
  it('a node below a Null dependent is evaluated with the cells that waited for it', async () => {
    const state = newState();
    await commit(state, (b) => {
      b.toRoot().apply(Node, { value: 1 }, { ref: 'source' });
      b.toRoot().apply(Gate, {}, { ref: 'gate', dependsOn: ['source'] });
    });
    // two edits: 'below' gets no cell
    await commit(state, (b) => {
      b.to('gate').apply(Node, { value: 1 }, { ref: 'below' });
      b.toRoot().apply(Node, { value: 1 }, { ref: 'other' });
    });
    await commit(state, (b) => {
      b.toRoot().apply(Node, { value: 1 }, { ref: 'user', dependsOn: ['below'] });
    });
    expect(state.cells.get('user')!.status).toBe('error');

    // 'gate' is a dependent of 'source' and gets an object; 'below' is made and 'user' evaluated in the same update
    await commit(state, (b) => {
      b.to('source').update({ value: 5 });
    });
    expect(valueOf(state, 'gate')).toBe(5);
    expect(valueOf(state, 'below')).toBe(6);
    expect(state.cells.get('user')!.status).toBe('ok');
    expect(valueOf(state, 'user')).toBe(7);
  });

  it('a dependent of a deleted dependency is evaluated again, unresolved until a cell of that ref is back', async () => {
    const state = newState();
    await commit(state, (b) => {
      b.toRoot().apply(Node, { value: 2 }, { ref: 'r' });
      b.toRoot().apply(Uses, { value: 1, source: 'r' }, { ref: 'user' });
    });
    expect(valueOf(state, 'user')).toBe(3);
    const user = state.cells.get('user')!;

    // a dependency derived from params is not deleted with its source; its object was made from it
    await commit(state, (b) => b.delete('r'));
    expect(state.cells.has('user')).toBe(true);
    expect(user.dependencies.dependsOn).toEqual([]);
    expect(user.status).toBe('error');
    expect(user.obj).toBeUndefined();

    await commit(state, (b) => b.toRoot().apply(Node, { value: 4 }, { ref: 'r' }));
    expect(user.status).toBe('ok');
    expect(valueOf(state, 'user')).toBe(5);
    expect(state.cells.get('r')!.dependencies.dependentBy).toEqual([user]);
  });
  it('forgets the waits of cells that no longer depend or are deleted', async () => {
    const state = newState();
    await commit(state, (b) => {
      b.toRoot().apply(Node, { value: 0 }, { ref: 'null' });
    });
    await commit(state, (b) => {
      b.to('null').apply(Node, { value: 1 }, { ref: 'below' });
      b.toRoot().apply(Node, { value: 1 }, { ref: 'other' });
    });
    await commit(state, (b) => {
      b.toRoot().apply(Uses, { value: 1, source: 'below' }, { ref: 'a' });
      b.toRoot().apply(Uses, { value: 1, source: 'below' }, { ref: 'b' });
    });
    expect(waiters(state).byCell.size).toBe(2);

    await commit(state, (b) => {
      b.to('a').update({ value: 1, source: '' });
      b.delete('b');
    });
    expect(state.cells.get('a')!.status).toBe('ok');
    expect(waiters(state).byCell.size).toBe(0);
  });
  it('a dependency cleared below a new Null object fails its dependent, not the update', async () => {
    const state = newState();
    await commit(state, (b) => {
      b.toRoot().apply(Node, { value: 1 }, { ref: 'parent' }).apply(Node, { value: 1 }, { ref: 'source' });
      b.toRoot().apply(Uses, { value: 1, source: 'source' }, { ref: 'user' });
    });
    expect(valueOf(state, 'user')).toBe(3);

    // 'source' has no object this pass; 'user' is not deferred for it
    await commit(state, (b) => b.to('parent').update({ value: 0 }));
    expect(state.cells.get('parent')!.obj).toBe(StateObject.Null);
    expect(state.cells.get('source')!.obj).toBeUndefined();
    expect(state.cells.get('user')!.status).toBe('error');
    expect(state.cells.get('user')!.errorText).toBe('Error: Unresolved dependency.');
  });
  it('a dependency evaluated below a Null object without one fails its dependent, not the update', async () => {
    const state = newState();
    await commit(state, (b) => {
      b.toRoot().apply(Node, { value: 0 }, { ref: 'null' });
    });
    // one edit: evaluated below the Null object, with no object
    await commit(state, (b) => {
      b.to('null').apply(Node, { value: 1 }, { ref: 'source' });
    });
    await commit(state, (b) => {
      b.toRoot().apply(Uses, { value: 1, source: 'source' }, { ref: 'user' });
    });
    expect(state.cells.get('user')!.status).toBe('error');

    // 'source' is evaluated again (status 'ok', still no object); 'user' is a dependent and not deferred for it
    await commit(state, (b) => b.to('source').update({ value: 2 }));
    expect(state.cells.get('source')!.status).toBe('ok');
    expect(state.cells.get('user')!.status).toBe('error');
  });
});
