/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { describe, expect, it } from 'vitest';
import { State, StateObject, StateTransform, StateTransformer, StateTree } from '../index.js';
import { Task } from '@molstar/core/task';
import { ParamDefinition as PD } from '@molstar/core/util/param-definition';

interface TypeInfo {
  name: string;
  typeClass: 'Root' | 'Data';
}
const Create = StateObject.factory<TypeInfo>();

class Root extends Create({ name: 'Root', typeClass: 'Root' }) {}
class Leaf extends Create<{ value: number }>({ name: 'Leaf', typeClass: 'Data' }) {}

const NS = 'state-cell-state-spec';

const CreateLeaf = StateTransformer.create<Root, Leaf, { value: number }>(NS, {
  name: 'create-leaf',
  from: [Root],
  to: [Leaf],
  display: { name: 'Create Leaf' },
  params: () => ({ value: PD.Numeric(0) }) as any,
  apply({ params }) {
    return new Leaf({ value: params.value });
  },
  update() {
    return StateTransformer.UpdateResult.Recreate;
  },
});

function newState() {
  return State.create(new Root({}), { runTask: <T>(t: Task<T>) => t.run() });
}

async function withLeaves(...refs: string[]) {
  const state = newState();
  const b = state.build();
  for (const ref of refs) b.toRoot().apply(CreateLeaf, { value: 0 }, { ref });
  await state.runTask(state.updateTree(b));
  return state;
}

function treeHidden(tree: StateTree, ref: string) {
  return !!tree.transforms.get(ref).state.isHidden;
}

function snapshotHidden(state: State, ref: string) {
  const json = state.getSnapshot().tree.transforms.find((t) => t.ref === ref)!;
  return !!json.state?.isHidden;
}

// tree.md §2, §4: flags written by updateCellState reach state.tree whatever builders are outstanding
describe('state flags and outstanding builders', () => {
  it('a flag set after build() reaches the tree', async () => {
    const state = await withLeaves('a', 'b');
    state.updateCellState('a', { isHidden: true });
    state.build(); // never committed
    state.updateCellState('b', { isHidden: true });

    expect(state.cells.get('b')!.state.isHidden).toBe(true);
    expect(treeHidden(state.tree, 'a')).toBe(true);
    expect(treeHidden(state.tree, 'b')).toBe(true);
    expect(snapshotHidden(state, 'b')).toBe(true);
  });

  it('a flag set after a transaction snapshot reaches the tree', async () => {
    const state = await withLeaves('a', 'b');
    state.updateCellState('a', { isHidden: true });
    await state.runTask(
      state.transaction(() => {
        state.updateCellState('b', { isHidden: true });
      }),
    );
    expect(treeHidden(state.tree, 'b')).toBe(true);
  });

  it('a builder made before a flag write keeps the flag on commit', async () => {
    const state = await withLeaves('a', 'b');
    const b = state.build();
    b.to('a').update({ value: 1 });
    state.updateCellState('b', { isHidden: true });
    await state.runTask(state.updateTree(b));

    expect(state.cells.get('a')!.obj!.data.value).toBe(1);
    expect(state.cells.get('b')!.state.isHidden).toBe(true);
    expect(treeHidden(state.tree, 'b')).toBe(true);
    expect(snapshotHidden(state, 'b')).toBe(true);
  });

  it('a builder with several edits keeps the flag on commit', async () => {
    const state = await withLeaves('a', 'b');
    const b = state.build();
    b.to('a').update({ value: 1 });
    b.toRoot().apply(CreateLeaf, { value: 2 }, { ref: 'c' });
    state.updateCellState('b', { isHidden: true });
    await state.runTask(state.updateTree(b));
    expect(treeHidden(state.tree, 'b')).toBe(true);
  });

  it('a builder rebased on a newer tree keeps the flag', async () => {
    const state = await withLeaves('a', 'b');
    const b = state.build();
    b.to('a').update({ value: 1 });
    const other = state.build();
    other.toRoot().apply(CreateLeaf, { value: 3 }, { ref: 'c' });
    await state.runTask(state.updateTree(other));
    state.updateCellState('b', { isHidden: true });
    await state.runTask(state.updateTree(b));
    expect(treeHidden(state.tree, 'b')).toBe(true);
    expect(state.cells.has('c')).toBe(true);
  });

  it('a flag survives a save and restore after a stale builder commit', async () => {
    const state = await withLeaves('a', 'b');
    const b = state.build();
    b.to('a').update({ value: 1 });
    state.updateCellState('b', { isHidden: true });
    await state.runTask(state.updateTree(b));

    const snapshot = state.getSnapshot();
    state.updateCellState('b', { isHidden: false });
    await state.runTask(state.setSnapshot(snapshot));
    expect(state.cells.get('b')!.state.isHidden).toBe(true);
  });

  it('undo restores the flags of the tree it saved', async () => {
    const state = await withLeaves('a', 'b');
    const b = state.build();
    b.to('a').update({ value: 1 });
    state.updateCellState('b', { isHidden: true });
    await state.runTask(state.updateTree(b));

    await state.runTask(state.updateTree(state.build().to('a').update({ value: 2 }), { canUndo: true }));
    await state.runTask(state.undo());
    expect(state.cells.get('a')!.obj!.data.value).toBe(1);
    expect(state.cells.get('b')!.state.isHidden).toBe(true);
  });

  it('a builder made before a flag write keeps its own state edits', async () => {
    const state = await withLeaves('a', 'b');
    const b = state.build();
    b.to('a').updateState({ isCollapsed: true });
    state.updateCellState('b', { isHidden: true });
    const tree = b.getTree();
    expect(tree.transforms.get('a').state.isCollapsed).toBe(true);
    await state.runTask(state.updateTree(b));
    expect(state.tree.transforms.get('a').state.isCollapsed).toBe(true);
    expect(treeHidden(state.tree, 'b')).toBe(true);
  });

  it('a builder keeps its own flag on a node whose other flag was set after build()', async () => {
    const state = await withLeaves('a');
    const b = state.build();
    b.to('a').updateState({ isCollapsed: true });
    state.updateCellState('a', { isHidden: true });
    await state.runTask(state.updateTree(b));
    expect(state.tree.transforms.get('a').state).toEqual({ isCollapsed: true, isHidden: true });
  });

  it('carries only the flags written after build(), not the rest of the cell', async () => {
    const state = await withLeaves('a');
    // a committed builder sets the tree's flag only (cells take flags from builders by updateCellState)
    const x = state.build();
    x.to('a').updateState({ isCollapsed: true });
    await state.runTask(state.updateTree(x));
    expect(state.cells.get('a')!.state.isCollapsed).toBeUndefined();

    const y = state.build();
    state.updateCellState('a', { isHidden: true });
    await state.runTask(state.updateTree(y));
    expect(state.tree.transforms.get('a').state).toEqual({ isCollapsed: true, isHidden: true });
  });
});

describe('TransientTree.asImmutable', () => {
  it('a builder can be edited after getTree()', () => {
    const state = newState();
    const b = state.build();
    b.toRoot().apply(CreateLeaf, { value: 0 }, { ref: 'a' });
    b.getTree();
    b.toRoot().apply(CreateLeaf, { value: 1 }, { ref: 'b' });
    const tree = b.getTree();
    expect(tree.transforms.has('a')).toBe(true);
    expect(tree.transforms.has('b')).toBe(true);
    expect(tree.children.get(StateTransform.RootRef).toArray()).toEqual(['a', 'b']);
  });

  it('a snapshot does not see later flag writes', async () => {
    const state = await withLeaves('a');
    state.updateCellState('a', { isHidden: true });
    const snapshot = state.build().getTree();
    state.updateCellState('a', { isCollapsed: true });
    expect(snapshot.transforms.get('a').state.isHidden).toBe(true);
    expect(snapshot.transforms.get('a').state.isCollapsed).toBeUndefined();
    expect(state.tree.transforms.get('a').state.isCollapsed).toBe(true);
  });
});
