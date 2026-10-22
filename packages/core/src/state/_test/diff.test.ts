/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { beforeEach, describe, expect, it, vi } from 'vitest';
import { State, StateObject, StateTransformer } from '../index.js';
import { findDeletes, findEditedDeletes, findEditedRoots, findUpdateRoots } from '../state/diff.js';
import { Task } from '@molstar/core/task';
import { ParamDefinition as PD } from '@molstar/core/util/param-definition';

// A builder on the current tree that moved no node finds its deletes and update roots from its edits (reconcile.md §2); they must be those of the whole-tree diff, in the same order.

interface TypeInfo {
  name: string;
  typeClass: 'Root' | 'Data';
}
const Create = StateObject.factory<TypeInfo>();
class Root extends Create({ name: 'Root', typeClass: 'Root' }) {}
class Leaf extends Create<{ value: number }>({ name: 'Leaf', typeClass: 'Data' }) {}

/** value 0: a Null object; negative: an error */
const Node = StateTransformer.create<Root | Leaf, Leaf, { value: number }>('update-roots-test', {
  name: 'node',
  from: [Root, Leaf],
  to: [Leaf],
  display: { name: 'Node' },
  params: () => ({ value: PD.Numeric(1) }),
  apply({ a, params }) {
    if (params.value === 0) return StateObject.Null;
    if (params.value < 0) throw new Error('negative');
    return new Leaf({ value: ((a.data as { value?: number }).value ?? 0) + params.value });
  },
});

function newState() {
  return State.create(new Root({}), { runTask: <T>(t: Task<T>) => t.run() });
}

/** a small deterministic generator (mulberry32) */
function random(seed: number) {
  return () => {
    seed = (seed + 0x6d2b79f5) | 0;
    let t = Math.imul(seed ^ (seed >>> 15), 1 | seed);
    t = (t + Math.imul(t ^ (t >>> 7), 61 | t)) ^ t;
    return ((t ^ (t >>> 14)) >>> 0) / 4294967296;
  };
}

describe('the diff from the edits of a builder', () => {
  beforeEach(() => {
    vi.spyOn(console, 'error').mockImplementation(() => {});
  });

  it('is the whole-tree diff, in the same order, over random edits', async () => {
    let compared = 0;
    let withRemoval = 0;
    let withRemovals = 0;
    for (let seed = 1; seed <= 40; seed++) {
      const rnd = random(seed);
      const pick = <T>(xs: readonly T[]) => xs[Math.floor(rnd() * xs.length)];
      const state = newState();
      // explicit refs: a seed picks the same nodes on every run (generated refs and their collisions: transient.test.ts)
      let n = 0;
      // values: mostly ok, some Null objects (0) and errors (-1)
      const value = () => pick([1, 2, 3, 4, 0, -1]);
      for (let step = 0; step < 60; step++) {
        const b = state.build();
        const refs = () => Array.from(b.currentTree.transforms.keys());
        const count = 1 + Math.floor(rnd() * 5);
        for (let e = 0; e < count; e++) {
          const nonRoot = refs().filter((r) => r !== state.tree.root.ref);
          const op = rnd();
          if (op < 0.45 || nonRoot.length === 0) b.to(pick(refs())).apply(Node, { value: value() }, { ref: `n${n++}` });
          else if (op < 0.65) b.to(pick(nonRoot)).update({ value: value() });
          else if (op < 0.75) b.to(pick(nonRoot)).updateState({ isHidden: rnd() < 0.5 });
          else if (op < 0.83) b.to(pick(nonRoot)).tag(`t${Math.floor(rnd() * 3)}`);
          else if (op < 0.88) b.currentTree.setTags(pick(nonRoot), [`d${Math.floor(rnd() * 3)}`]);
          else if (op < 0.92)
            b.to(pick(nonRoot)).update({ value: 1 }); // often unchanged
          else if (op < 0.97) b.delete(pick(nonRoot));
          else {
            // a node and one of its ancestors, or two unrelated ones, in either order
            const r = pick(nonRoot);
            const t = b.currentTree.transforms.get(r);
            b.delete(rnd() < 0.5 ? r : pick(nonRoot));
            if (b.currentTree.transforms.has(t.parent)) b.delete(t.parent);
          }
        }
        const tree = b.getTree();
        const info = b.editInfo;
        // the cases the reconciler takes the edits for: not the fast track, nothing moved, cells in line with the tree
        const inSync = (state as unknown as { cellsInSync: boolean }).cellsInSync;
        const edits = info.edits!;
        if (inSync && !edits.moved && info.count !== 1) {
          const cells = state.cells;
          expect(findEditedDeletes(state.tree, tree, cells, edits.removed)).toEqual(
            findDeletes(state.tree, tree, cells),
          );
          expect(findEditedRoots(cells, tree, edits.refs)).toEqual(findUpdateRoots(cells, tree));
          compared++;
          if (edits.removed.length) withRemoval++;
          if (edits.removed.length > 1) withRemovals++;
        }
        await state.runTask(state.updateTree(b, { doNotLogTiming: true }));
      }
    }
    expect(compared).toBeGreaterThan(300);
    expect(withRemoval).toBeGreaterThan(20);
    expect(withRemovals).toBeGreaterThan(20);
  });

  it('finds the cells an update that threw left behind, as the whole-tree diff does', async () => {
    const state = newState();
    const b0 = state.build();
    b0.toRoot().apply(Node, { value: 1 }, { ref: 'a' });
    b0.toRoot().apply(Node, { value: 1 }, { ref: 'b' });
    b0.toRoot().apply(Node, { value: 1 }, { ref: 'x' });
    await state.runTask(state.updateTree(b0));

    // a dependency on a missing ref throws while the cells are made: 'x' has its new params in the tree, its cell the
    // old transform
    const failing = state.build();
    failing.to('x').update({ value: 5 });
    failing.toRoot().apply(Node, { value: 1 }, { ref: 'bad', dependsOn: ['missing'] });
    await expect(state.runTask(state.updateTree(failing))).rejects.toThrow(/non-existent transform/);
    expect(state.cells.get('x')!.transform.version).not.toBe(state.tree.transforms.get('x').version);

    // two edits elsewhere: not the fast track; 'x' is evaluated again
    const next = state.build();
    next.to('a').update({ value: 2 });
    next.to('b').update({ value: 3 });
    await state.runTask(state.updateTree(next));
    expect((state.cells.get('x')!.obj as Leaf).data.value).toBe(5);
  });

  it('finds the cells a fast track with an edit outside its node left behind', async () => {
    const state = newState();
    const b0 = state.build();
    b0.toRoot().apply(Node, { value: 1 }, { ref: 'a' });
    b0.toRoot().apply(Node, { value: 1 }, { ref: 'b' });
    b0.toRoot().apply(Node, { value: 1 }, { ref: 'x' });
    await state.runTask(state.updateTree(b0));

    // one counted edit takes the fast track; the tag written to the builder's tree directly is not evaluated
    const fast = state.build();
    fast.to('a').update({ value: 2 });
    fast.currentTree.setTags('x', ['tagged']);
    await state.runTask(state.updateTree(fast));
    expect(state.cells.get('x')!.transform.tags).toBeUndefined();

    const next = state.build();
    next.to('a').update({ value: 3 });
    next.to('b').update({ value: 3 });
    await state.runTask(state.updateTree(next));
    expect(state.cells.get('x')!.transform.tags).toEqual(['tagged']);
  });

  it('gives the same cells and events as the whole-tree diff, through updateTree', async () => {
    let fromEdits = 0;
    let withRemovals = 0;
    for (let seed = 1; seed <= 30; seed++) {
      const rnd = random(1000 + seed);
      const pick = <T>(xs: readonly T[]) => xs[Math.floor(rnd() * xs.length)];
      // `twin` takes the whole trees every time: cellsInSync is cleared before each of its commits
      const [state, twin] = [newState(), newState()];
      const log = new Map<State, string[]>([
        [state, []],
        [twin, []],
      ]);
      for (const s of [state, twin]) {
        const l = log.get(s)!;
        s.events.object.created.subscribe((e) => l.push(`+${e.ref}`));
        s.events.object.updated.subscribe((e) => l.push(`~${e.ref}:${e.action}`));
        s.events.object.removed.subscribe((e) => l.push(`-${e.ref}`));
        s.events.cell.created.subscribe((e) => l.push(`c+${e.ref}`));
        s.events.cell.removed.subscribe((e) => l.push(`c-${e.ref}`));
      }
      let n = 0;
      for (let step = 0; step < 40; step++) {
        // the same edits, decided once, on both states (explicit refs: the trees have the same refs and order)
        const ops: ((b: ReturnType<State['build']>) => void)[] = [];
        const refs = Array.from(state.tree.transforms.keys());
        const added: string[] = [];
        const count = 1 + Math.floor(rnd() * 4);
        for (let e = 0; e < count; e++) {
          const all = [...refs, ...added];
          const nonRoot = all.filter((r) => r !== state.tree.root.ref);
          const op = rnd();
          if (op < 0.4 || nonRoot.length === 0) {
            const parent = pick(all);
            const ref = `n${n++}`;
            const value = pick([1, 2, 0, -1]);
            const dependsOn =
              rnd() < 0.15 && refs.length > 1 ? [pick(refs.filter((r) => r !== state.tree.root.ref) || [])] : undefined;
            added.push(ref);
            ops.push((b) => {
              if (b.currentTree.transforms.has(parent))
                b.to(parent).apply(
                  Node,
                  { value },
                  { ref, dependsOn: dependsOn?.filter((d) => d && b.currentTree.transforms.has(d)) },
                );
            });
          } else if (op < 0.7) {
            const r = pick(nonRoot);
            const value = pick([1, 2, 3, 0, -1]);
            ops.push((b) => void (b.currentTree.transforms.has(r) && b.to(r).update({ value })));
          } else if (op < 0.8) {
            const r = pick(nonRoot);
            const isHidden = rnd() < 0.5;
            ops.push((b) => void (b.currentTree.transforms.has(r) && b.to(r).updateState({ isHidden })));
          } else {
            const r = pick(nonRoot);
            ops.push((b) => void (b.currentTree.transforms.has(r) && b.delete(r)));
          }
        }
        // sometimes a builder made before another commit (rebased when it commits)
        const rebase = rnd() < 0.15;
        for (const s of [state, twin]) {
          const b = s.build();
          for (const op of ops) op(b);
          if (rebase) {
            const other = s.build();
            other.toRoot().apply(Node, { value: 1 }, { ref: `r${step}` });
            if (s === twin) (s as unknown as { cellsInSync: boolean }).cellsInSync = false;
            await s.runTask(s.updateTree(other, { doNotLogTiming: true })).catch(() => {});
          }
          if (s === twin) (s as unknown as { cellsInSync: boolean }).cellsInSync = false;
          else if (!rebase && b.editInfo.count !== 1 && (s as unknown as { cellsInSync: boolean }).cellsInSync) {
            fromEdits++;
            // several removals, or one that removed its dependents too
            if (b.editInfo.edits!.removed.length > 1) withRemovals++;
          }
          await s.runTask(s.updateTree(b, { doNotLogTiming: true })).catch(() => {});
        }
        const cells = (s: State) =>
          Array.from(s.cells.values())
            .map((c) => {
              const t = s.tree.transforms.get(c.transform.ref);
              const value = c.obj === StateObject.Null ? 'null' : (c.obj?.data as { value?: number })?.value;
              return `${c.transform.ref}:${c.status}:${value}:${t ? t.version === c.transform.version : 'gone'}`;
            })
            .sort();
        expect(cells(state)).toEqual(cells(twin));
        expect(log.get(state)).toEqual(log.get(twin));
      }
    }
    expect(fromEdits).toBeGreaterThan(200);
    expect(withRemovals).toBeGreaterThan(20);
  });
});
