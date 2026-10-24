/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { describe, expect, it } from 'vitest';
import { State, StateObject, StateTransformer } from '../index.js';
import type { Task } from '@molstar/core/task';
import { ParamDefinition as PD } from '@molstar/core/util/param-definition';

interface TypeInfo {
  name: string;
  typeClass: 'Root' | 'Data';
}
const Create = StateObject.factory<TypeInfo>();
class Root extends Create({ name: 'Root', typeClass: 'Root' }) {}
class Leaf extends Create<{ value: number }>({ name: 'Leaf', typeClass: 'Data' }) {}

const Node = StateTransformer.create<Root, Leaf, { value: number }>('state-builder-tagged-test', {
  name: 'node',
  from: [Root],
  to: [Leaf],
  display: { name: 'Node' },
  params: () => ({ value: PD.Numeric(1) }),
  apply({ params }) {
    return new Leaf({ value: params.value });
  },
});

describe('applyOrUpdateTagged', () => {
  it('adds the extra tags when the params change too (perf-audit-2 03/B12)', async () => {
    const state = State.create(new Root({}), { runTask: <T>(t: Task<T>) => t.run() });
    await state.build().toRoot().applyOrUpdateTagged('a', Node, { value: 1 }).commit();

    const b = state.build();
    const ref = b.toRoot().applyOrUpdateTagged('a', Node, { value: 2 }, { tags: ['extra'] }).ref;
    await b.commit();

    const cell = state.cells.get(ref)!;
    expect(cell.transform.params).toEqual({ value: 2 });
    expect(cell.transform.tags).toEqual(['a', 'extra']);
    expect(cell.obj?.data.value).toBe(2);
    expect(cell.obj?.tags).toEqual(['a', 'extra']);
  });
});

describe('rebased builder', () => {
  it('replays tags and dependsOn onto the newer tree (perf-audit-2 03/B3)', async () => {
    const state = State.create(new Root({}), { runTask: <T>(t: Task<T>) => t.run() });
    const b0 = state.build();
    b0.toRoot().apply(Node, { value: 1 }, { ref: 'a', tags: 'a' });
    b0.toRoot().apply(Node, { value: 1 }, { ref: 'b' });
    await b0.commit();

    // made on the current tree, committed after another commit: rebased
    const b = state.build();
    b.to('a').tag('extra');
    b.to('b').dependsOn('a');
    await state.build().toRoot().apply(Node, { value: 3 }, { ref: 'c' }).commit();
    await b.commit();

    expect(state.tree.transforms.get('a').tags).toEqual(['a', 'extra']);
    expect(state.tree.transforms.get('b').dependsOn).toEqual(['a']);
    expect(state.tree.dependencies.get('a')?.toArray()).toEqual(['b']);
    expect(state.cells.get('a')!.obj?.tags).toEqual(['a', 'extra']);
    expect(state.cells.has('c')).toBe(true);
  });

  it("adds the replayed tags and dependsOn to the newer tree's", async () => {
    const state = State.create(new Root({}), { runTask: <T>(t: Task<T>) => t.run() });
    const b0 = state.build();
    b0.toRoot().apply(Node, { value: 1 }, { ref: 'a', tags: 'a' });
    b0.toRoot().apply(Node, { value: 1 }, { ref: 'b' });
    b0.toRoot().apply(Node, { value: 1 }, { ref: 'c' });
    await b0.commit();

    const b = state.build();
    b.to('a').tag('x');
    b.to('c').dependsOn('a');
    // committed first: another tag and dependency
    const other = state.build();
    other.to('a').tag('z');
    other.to('c').dependsOn('b');
    await other.commit();
    await b.commit();

    expect(state.tree.transforms.get('a').tags).toEqual(['a', 'x', 'z']);
    expect(state.tree.transforms.get('c').dependsOn).toEqual(['a', 'b']);
  });

  it('replays the extra tags of applyOrUpdateTagged with unchanged params', async () => {
    const state = State.create(new Root({}), { runTask: <T>(t: Task<T>) => t.run() });
    await state.build().toRoot().applyOrUpdateTagged('a', Node, { value: 1 }).commit();

    const b = state.build();
    const ref = b.toRoot().applyOrUpdateTagged('a', Node, { value: 1 }, { tags: ['extra'] }).ref;
    await state.build().toRoot().apply(Node, { value: 3 }).commit();
    await b.commit();
    expect(state.tree.transforms.get(ref).tags).toEqual(['a', 'extra']);
  });
});
