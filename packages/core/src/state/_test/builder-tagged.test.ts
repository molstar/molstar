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
