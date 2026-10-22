/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Frames of an animation applied one after another with `setSnapshot`, as a playing snapshot transition does: each is
// a serialized tree, and consecutive frames share the serialized transforms that did not change (MVS builds them so).

import { describe, expect, it, vi } from 'vitest';
import { Task } from '@molstar/core/task';
import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import { State, StateObject, StateTransform, StateTransformer, type StateTree } from '../index.js';

const Create = StateObject.factory<{ name: string; typeClass: 'Root' | 'Data' }>();
class Root extends Create({ name: 'Root', typeClass: 'Root' }) {}
class Value extends Create<{ value: number }>({ name: 'Value', typeClass: 'Data' }) {}

let counter = 0;
/** A transformer from `from` with a numeric param; `params` and `update` are spies */
function valueTransformer(from: typeof Root | typeof Value) {
  const params = vi.fn((_a?: unknown, _ctx?: unknown) => ({ value: PD.Numeric(0), extra: PD.Boolean(false) }));
  const update = vi.fn(({ b, newParams }: any) => {
    b.data.value = newParams.value;
    return StateTransformer.UpdateResult.Updated;
  });
  const transformer = StateTransformer.create<any, Value, { value: number; extra: boolean }>('snapshot-frames', {
    name: `value-${counter++}`,
    from: [from],
    to: [Value],
    display: { name: 'Value' },
    params: params as any,
    apply: ({ params }) => new Value({ value: params.value }),
    update,
  });
  return { transformer, params, update };
}

function serialized(
  transformer: StateTransformer,
  ref: string,
  parent: string,
  value: number,
): StateTree.Serialized['transforms'][number] {
  return { transformer: transformer.id, ref, parent, params: { value }, version: `${ref}:${value}` } as any;
}

/** frames 0..n-1: `a` takes the frame's value, `b` below it keeps one serialized transform for all frames */
function frames(n: number) {
  const a = valueTransformer(Root);
  const b = valueTransformer(Value);
  const root = { transformer: StateTransformer.ROOT.id, ref: StateTransform.RootRef, parent: StateTransform.RootRef };
  const bNode = serialized(b.transformer, 'b', 'a', 7);
  const snapshots: State.Snapshot[] = [];
  for (let i = 0; i < n; i++) {
    snapshots.push({
      tree: { transforms: [root as any, serialized(a.transformer, 'a', StateTransform.RootRef, i), bNode] },
    });
  }
  return { a, b, snapshots };
}

const newState = () => State.create(new Root({}), { runTask: <T>(t: Task<T>) => t.run() });

describe('frames applied with setSnapshot', () => {
  it('evaluate the changed node and its subtree', async () => {
    const state = newState();
    const { b, snapshots } = frames(4);
    for (const s of snapshots) await state.setSnapshot(s).run();
    expect(state.cells.get('a')!.obj!.data.value).toBe(3);
    // b's transform is the same in every frame, but its parent changed: it is evaluated each time
    expect(b.update).toHaveBeenCalledTimes(3);
  });

  it('build the default param definition of a transformer once (S9)', async () => {
    const state = newState();
    const { a, b, snapshots } = frames(4);
    for (const s of snapshots) await state.setSnapshot(s).run();
    // `params(undefined, undefined)`: what dependencies are found from (StateTransform.getEffectiveDependsOn)
    const withoutSource = (p: typeof a.params) => p.mock.calls.filter((c) => c[0] === undefined).length;
    expect(withoutSource(a.params)).toBe(1);
    expect(withoutSource(b.params)).toBe(1);
  });
});
