/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Map as ImmutableMap } from 'immutable';
import { describe, expect, it, vi } from 'vitest';
import { State, StateObject, StateTransform, StateTransformer, StateTree } from '../index.js';
import { Task } from '@molstar/core/task';
import { UUID } from '@molstar/core/util/uuid';
import { ParamDefinition as PD } from '@molstar/core/util/param-definition';

// A transient tree's edits never reach the tree it was made on (tree.md §2), whatever it removed before them.

interface TypeInfo {
  name: string;
  typeClass: 'Root' | 'Data';
}
const Create = StateObject.factory<TypeInfo>();
class Root extends Create({ name: 'Root', typeClass: 'Root' }) {}
class Leaf extends Create<{ value: number }>({ name: 'Leaf', typeClass: 'Data' }) {}

/** value 0: a Null object; negative: an error */
const Node = StateTransformer.create<Root | Leaf, Leaf, { value: number }>('state-transient-test', {
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

/**
 * Refs in pairs with equal Immutable hashes ('Aa' and 'BB' hash alike): with more than 8 keys a map is a hashed trie,
 * and the two share a collision node. Generated refs (UUIDs that differ in the last characters) collide this way often.
 */
const colliding = (n: number) => `n${n >> 1}${n & 1 ? 'Aa' : 'BB'}`;

const Alphabet = 'ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789-_';
/**
 * `UUID.create22`'s layout (a 14-character prefix, then an 8-character base64 counter) with a fixed prefix: the
 * collisions of generated refs and versions, the same in every run, so a failing seed reproduces.
 */
function generatedIds() {
  let n = 0;
  const spy = vi.spyOn(UUID, 'create22').mockImplementation(() => {
    let s = '';
    for (let i = 0, c = n++; i < 8; i++, c = Math.floor(c / 64)) s = Alphabet[c % 64] + s;
    return `transientTests${s}` as UUID;
  });
  return { reset: () => void (n = 0), restore: () => spy.mockRestore() };
}

// The workaround in `TransientTree` (`renewOwner`, tree.md §2) is for this. When it fails, Immutable is fixed: see the
// core TODO.
it('Immutable writes a set after a delete of a colliding key into the frozen map', () => {
  const entries: [string, number][] = [
    ['Aa', 1],
    ['BB', 2],
  ];
  for (let i = 0; i < 10; i++) entries.push([`k${i}`, i]);
  const frozen = ImmutableMap(entries);
  const mutable = frozen.asMutable();
  mutable.delete('Aa');
  mutable.set('BB', 3);
  expect(frozen.get('BB')).toBe(3);
});

/** 10 nodes under the root, then the pair 'pAa' and 'pBB'; each with a dependent */
async function withPair() {
  const state = newState();
  const b = state.build();
  for (const ref of [...Array.from({ length: 10 }, (_, i) => `k${i}`), 'pAa', 'pBB']) {
    b.toRoot().apply(Node, { value: 1 }, { ref });
    b.toRoot().apply(Node, { value: 1 }, { ref: `${ref}-dep`, dependsOn: [ref] });
  }
  await state.runTask(state.updateTree(b));
  return state;
}

/** every node, its parent, version and children, and the dependencies: what a commit compares */
function shape(tree: StateTree) {
  const nodes: string[] = [];
  tree.transforms.forEach((t: StateTransform) => {
    const children = tree.children.get(t.ref);
    nodes.push(`${t.ref}<${t.parent}:${t.version}:[${children ? children.toArray().join(',') : 'none'}]`);
  });
  tree.children.forEach((_: unknown, ref: string) => {
    if (!tree.transforms.has(ref)) nodes.push(`children of a missing ${ref}`);
  });
  tree.dependencies.forEach((deps: StateTree.ChildSet, ref: string) => {
    nodes.push(`${ref}->[${deps.toArray().join(',')}]`);
  });
  return nodes.sort();
}

/** what breaks: each child set lists exactly the nodes whose parent it is, and each child has a transform */
function inconsistencies(tree: StateTree) {
  const found: string[] = [];
  tree.transforms.forEach((t: StateTransform) => {
    const children = tree.children.get(t.ref);
    if (!children) return void found.push(`no children of ${t.ref}`);
    children.forEach((c) => {
      if (tree.transforms.get(c!)?.parent !== t.ref) found.push(`child ${c} of ${t.ref}`);
    });
    if (t.ref !== t.parent && !tree.children.get(t.parent)?.has(t.ref)) found.push(`${t.ref} not in its parent`);
  });
  return found;
}

describe('a transient tree after a removal', () => {
  it('does not add a child to the tree it was made on', async () => {
    const state = await withPair();
    const before = shape(state.tree);
    const b = state.build();
    b.delete('pAa');
    b.to('pBB').apply(Node, { value: 1 }, { ref: 'c' });
    expect(state.tree.children.get('pBB').has('c')).toBe(false);
    expect(shape(state.tree)).toEqual(before);
    expect(inconsistencies(b.getTree())).toEqual([]);
  });

  it('does not change a transform of the tree it was made on', async () => {
    const state = await withPair();
    const before = state.tree.transforms.get('pBB');
    const b = state.build();
    b.delete('pAa');
    b.to('pBB').update({ value: 2 });
    expect(state.tree.transforms.get('pBB')).toBe(before);
    expect(b.currentTree.transforms.get('pBB').params.value).toBe(2);
  });

  it('does not add a dependent to the tree it was made on', async () => {
    const state = await withPair();
    const before = shape(state.tree);
    const b = state.build();
    b.delete('pAa');
    b.toRoot().apply(Node, { value: 1 }, { ref: 'd', dependsOn: ['pBB'] });
    expect(state.tree.dependencies.get('pBB').has('d')).toBe(false);
    expect(shape(state.tree)).toEqual(before);
    expect(b.getTree().dependencies.get('pBB').toArray()).toEqual(['pBB-dep', 'd']);
  });

  it('leaves the state committable while that builder is outstanding', async () => {
    const state = await withPair();
    const outstanding = state.build();
    outstanding.delete('pAa');
    outstanding.to('pBB').apply(Node, { value: 1 }, { ref: 'c' });

    // two edits: the whole-tree diff, which walked the leaked child 'c' without a transform
    const other = state.build();
    other.to('k0').update({ value: 2 });
    other.to('k1').update({ value: 3 });
    await state.runTask(state.updateTree(other));
    expect(state.cells.has('c')).toBe(false);

    // rebased on the new tree (tree.md §4)
    await state.runTask(state.updateTree(outstanding));
    expect(state.cells.has('pAa')).toBe(false);
    expect(state.cells.get('c')!.obj).toBeInstanceOf(Leaf);
    expect(inconsistencies(state.tree)).toEqual([]);
  });

  it('keeps a child added to the root after the root was cleared', async () => {
    const state = await withPair();
    const b = state.build();
    b.delete(StateTransform.RootRef);
    b.toRoot().apply(Node, { value: 1 }, { ref: 'a' });
    b.toRoot().apply(Node, { value: 1 }, { ref: 'b' });
    const tree = b.getTree();
    expect(tree.children.get(StateTransform.RootRef).toArray()).toEqual(['a', 'b']);
    expect(inconsistencies(tree)).toEqual([]);
    expect(state.tree.children.get(StateTransform.RootRef).size).toBe(24);
  });
});

/** a small deterministic generator (mulberry32) */
function random(seed: number) {
  return () => {
    seed = (seed + 0x6d2b79f5) | 0;
    let t = Math.imul(seed ^ (seed >>> 15), 1 | seed);
    t = (t + Math.imul(t ^ (t >>> 7), 61 | t)) ^ t;
    return ((t ^ (t >>> 14)) >>> 0) / 4294967296;
  };
}

describe('random builders with removals', () => {
  it.each([
    ['colliding refs', 'colliding'],
    ['generated refs', 'generated'],
  ] as const)('never change the state tree before they commit (%s)', async (_, refs) => {
    vi.spyOn(console, 'error').mockImplementation(() => {});
    const ids = generatedIds();
    let removals = 0;
    try {
      for (let seed = 1; seed <= 100; seed++) {
        ids.reset();
        const rnd = random(seed);
        const pick = <T>(xs: readonly T[]) => xs[Math.floor(rnd() * xs.length)];
        const state = newState();
        let n = 0;
        const options = () => (refs === 'colliding' ? { ref: colliding(n++) } : undefined);
        const value = () => pick([1, 2, 3, 4, 0, -1]);
        for (let step = 0; step < 60; step++) {
          const committed = shape(state.tree).join();
          const b = state.build();
          const all = () => Array.from(b.currentTree.transforms.keys());
          const count = 1 + Math.floor(rnd() * 5);
          for (let e = 0; e < count; e++) {
            const nonRoot = all().filter((r) => r !== StateTransform.RootRef);
            const op = rnd();
            if (op < 0.45 || nonRoot.length === 0) b.to(pick(all())).apply(Node, { value: value() }, options());
            else if (op < 0.65) b.to(pick(nonRoot)).update({ value: value() });
            else if (op < 0.75) b.to(pick(nonRoot)).tag(`t${Math.floor(rnd() * 3)}`);
            else if (op < 0.85) {
              const deps = nonRoot.filter(() => rnd() < 0.3).slice(0, 2);
              b.to(pick(all())).apply(Node, { value: 1 }, { ...options(), dependsOn: deps });
            } else {
              b.delete(pick(nonRoot));
              removals++;
            }
          }
          const where = `seed ${seed}, step ${step}`;
          const unchanged = shape(state.tree).join() === committed;
          if (!unchanged) expect(shape(state.tree), where).toEqual(committed.split(','));
          const built = inconsistencies(b.getTree());
          if (built.length) expect(built, where).toEqual([]);
          await state.runTask(state.updateTree(b, { doNotLogTiming: true }));
          const after = inconsistencies(state.tree);
          if (after.length) expect(after, where).toEqual([]);
        }
      }
    } finally {
      ids.restore();
    }
    expect(removals).toBeGreaterThan(500);
  });
});
