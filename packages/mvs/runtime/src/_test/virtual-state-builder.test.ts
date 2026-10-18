/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { describe, expect, it } from 'vitest';
import { PluginContext } from '@molstar/plugin/context';
import { StateTree } from '@molstar/core/state';
import { Download } from '@molstar/plugin/state/transforms/data/fetch';
import { ParseCif } from '@molstar/plugin/state/formats/cif';
import {
  CustomStructureProperties,
  StructureFromModel,
  TransformStructureConformation,
} from '@molstar/plugin/state/transforms/structure/hierarchy';
import { StructureRepresentation3D } from '@molstar/plugin/state/transforms/structure/representation';
import { jsonHash } from '../helpers/json-hash.js';
import {
  jsonEqual,
  type LoadingActions,
  loadTreeVirtual,
  shareEqualParts,
  SubtreeMemo,
  UpdateTarget,
  VirtualStateBuilder,
} from '../load-generic.js';

const Transform = { name: 'components', params: { translation: [1, 2, 3], angle: 0, axis: [0, 0, 1] } } as any;

/** The same edits on a `StateBuilder` and on a `VirtualStateBuilder`. */
function edits(api: {
  apply(parent: string, tr: any, params: any, options: any): void;
  tag(ref: string, tags: string[]): void;
  dependsOn(ref: string, refs: string[]): void;
  updateState(ref: string, state: any): void;
  root: string;
}) {
  api.apply(api.root, Download, { url: 'https://example.org/1cbs.bcif', isBinary: true }, { ref: 'a' });
  api.apply('a', ParseCif, {}, { ref: 'b', tags: ['z', 'y'] });
  api.apply('b', StructureFromModel, undefined, { ref: 's' });
  api.apply('s', StructureRepresentation3D, { type: { name: 'cartoon', params: {} } }, { ref: 'c' });
  api.apply('s', TransformStructureConformation, { transform: Transform }, { ref: 'd' }); // between s and c
  api.apply('s', CustomStructureProperties, { properties: {} }, { ref: 'e' }); // a decorator too: between s and d
  api.tag('c', ['x', 'a']);
  api.tag('c', ['a']);
  api.dependsOn('e', ['a']);
  api.dependsOn('e', ['a']);
  api.updateState('c', { isHidden: true });
  api.apply('d', TransformStructureConformation, { transform: Transform }, { ref: 'f' }); // between d and its children
  api.apply('s', StructureRepresentation3D, { type: { name: 'spacefill', params: {} } }, { ref: 'g' }); // under f, the decorator root of s
  api.apply(api.root, Download, { url: 'https://example.org/2.bcif', isBinary: true }, { ref: 'h', tags: 'one' });
}

describe('VirtualStateBuilder', () => {
  it('serializes as StateTree.toJSON of a StateBuilder tree with hash versions', async () => {
    const plugin = new PluginContext({ registry: [], behaviors: [] });
    await plugin.init();

    const builder = plugin.build();
    edits({
      root: builder.toRoot().ref,
      apply: (parent, tr, params, options) => builder.to(parent).apply(tr, params, options),
      tag: (ref, tags) => builder.to(ref).tag(tags),
      dependsOn: (ref, refs) => builder.to(ref).dependsOn(refs),
      updateState: (ref, state) => builder.to(ref).updateState(state),
    });
    // (the same transforms; the versions of those added are `jsonHash` of their params instead of `hashMurmur128o`)
    const root = plugin.state.data.tree.root.ref;
    const expected = {
      transforms: StateTree.toJSON(builder.getTree({ useHashVersion: true })).transforms.map((t) =>
        t.ref === root ? t : { ...t, version: jsonHash(t.params ?? {}, new WeakMap()) },
      ),
    };

    const virtual = new VirtualStateBuilder(plugin.state.data, true);
    edits({
      root: virtual.rootRef,
      apply: (parent, tr, params, options) => virtual.apply(parent, tr, params, options),
      tag: (ref, tags) => virtual.tag(ref, tags),
      dependsOn: (ref, refs) => virtual.dependsOn(ref, refs),
      updateState: (ref, state) => virtual.updateState(ref, state),
    });
    const actual = virtual.serialize();

    expect(actual).toEqual(expected);
    expect(actual.transforms.map((t) => t.ref)).toEqual([
      plugin.state.data.tree.root.ref,
      'a',
      'b',
      's',
      'e',
      'd',
      'f',
      'c',
      'g',
      'h',
    ]);

    // reused: equal transforms of an earlier tree are those objects
    const again = new VirtualStateBuilder(plugin.state.data, true);
    edits({
      root: again.rootRef,
      apply: (parent, tr, params, options) => again.apply(parent, tr, params, options),
      tag: (ref, tags) => again.tag(ref, tags),
      dependsOn: (ref, refs) => again.dependsOn(ref, refs),
      updateState: (ref, state) => again.updateState(ref, state),
    });
    const reused = again.serialize([actual]);
    expect(reused).toEqual(expected);
    // (the root is not built here: serialized again, not looked up)
    expect(reused.transforms.map((t, i) => t === actual.transforms[i])).toEqual([
      false,
      true,
      true,
      true,
      true,
      true,
      true,
      true,
      true,
      true,
    ]);

    plugin.dispose();
  });
});

describe('loadTreeVirtual with a SubtreeMemo', () => {
  // a toy tree: every child of the root downloads its url, every grandchild parses
  const actions: LoadingActions<any, {}> = {
    root: (parent) => parent,
    file: (parent, node: any) => UpdateTarget.apply(parent, Download, { url: node.params.url, isBinary: true }),
    cif: (parent) => UpdateTarget.apply(parent, ParseCif, {}),
  } as any;
  const file = (url: string) => ({ kind: 'file', params: { url }, children: [{ kind: 'cif' }] });

  it('replays unchanged units with the same output, and loads again when a replay does not fit', async () => {
    const plugin = new PluginContext({ registry: [], behaviors: [] });
    await plugin.init();
    const load = (tree: any, memo?: SubtreeMemo) => {
      if (memo) memo.units = new Set(tree.children);
      return loadTreeVirtual(plugin, tree, actions, {}, { replaceExisting: true, subtreeMemo: memo }).data!.tree;
    };
    const memo = new SubtreeMemo(() => true);
    const a = file('https://example.org/a.bcif');
    const b = file('https://example.org/b.bcif');
    const stats0 = { ...SubtreeMemo.stats };

    const first = load({ kind: 'root', children: [a, b] }, memo);
    expect(first).toEqual(load({ kind: 'root', children: [a, b] }));

    // b changed: a is replayed, b loaded
    const b2 = file('https://example.org/b2.bcif');
    const second = load({ kind: 'root', children: [a, b2] }, memo);
    expect(second).toEqual(load({ kind: 'root', children: [a, b2] }));
    expect(second.transforms[1]).toBe(first.transforms[1]);
    expect(SubtreeMemo.stats.replayed - stats0.replayed).toBe(1);

    // a's child is a new object: a's own inputs are the same objects, and so are the child's (none): both replay their own
    // transforms, nothing is loaded
    const a2 = { ...a, children: [{ kind: 'cif' }] };
    const ownBefore = SubtreeMemo.stats.replayedOwn;
    const loadedBefore = SubtreeMemo.stats.loaded;
    const own = load({ kind: 'root', children: [a2, b2] }, memo);
    expect(own).toEqual(load({ kind: 'root', children: [a2, b2] }));
    expect(SubtreeMemo.stats.replayedOwn - ownBefore).toBe(2);
    expect(SubtreeMemo.stats.loaded - loadedBefore).toBe(0);
    expect(own.transforms[1]).toBe(first.transforms[1]);

    // a unit moved: at its position before, loaded again (no replay of another node's record)
    const c = file('https://example.org/c.bcif');
    const third = load({ kind: 'root', children: [c, a2, b2] }, memo);
    expect(third).toEqual(load({ kind: 'root', children: [c, a2, b2] }));

    // the same unit at the same position, but its ref counter is not the recorded one (the unit before it no longer
    // takes one): the replay does not fit, the tree is loaded again
    const fourth = load({ kind: 'root', children: [{ kind: 'note', children: [] }, a2, b2] }, memo);
    expect(fourth).toEqual(load({ kind: 'root', children: [{ kind: 'note', children: [] }, a2, b2] }));
    expect(SubtreeMemo.stats.fallbacks - stats0.fallbacks).toBe(1);
    expect(memo.fallbacks).toBe(1);

    plugin.dispose();
  });
});

describe('jsonEqual', () => {
  it('is JSON.stringify equality', () => {
    const cases: [unknown, unknown][] = [
      [
        { a: 1, b: [1, 2] },
        { a: 1, b: [1, 2] },
      ],
      [
        { a: 1, b: 2 },
        { b: 2, a: 1 },
      ],
      [{ a: 1, b: undefined }, { a: 1 }],
      [{ a: () => 1 }, {}],
      [[undefined], [null]],
      [
        [1, NaN],
        [1, null],
      ],
      [{ a: Infinity }, { a: null }],
      [-0, 0],
      [new Float32Array([1, 2]), { 0: 1, 1: 2 }],
      [{ d: new Date(0) }, { d: '1970-01-01T00:00:00.000Z' }],
      [{ a: [1, { b: 'x' }] }, { a: [1, { b: 'y' }] }],
      [1, '1'],
      [[], {}],
    ];
    for (const [a, b] of cases) {
      const same = JSON.stringify(a) === JSON.stringify(b);
      expect(jsonEqual(a, b), JSON.stringify([a, b])).toBe(same);
      // shareEqualParts decides equality in the same walk: the old value exactly when the JSON is equal
      expect(shareEqualParts(a, b) === b, JSON.stringify([a, b])).toBe(same);
      expect(JSON.stringify(shareEqualParts(a, b))).toBe(JSON.stringify(a));
    }
  });
});

describe('shareEqualParts', () => {
  it('keeps the JSON of the fresh value and takes the equal parts of the old one', () => {
    const old = {
      type: { name: 'cartoon', params: { alpha: 1, x: [1, 2] } },
      colorTheme: { name: 'uniform', params: { value: 3 } },
    };
    const fresh = {
      type: { name: 'cartoon', params: { alpha: 0.5, x: [1, 2] } },
      colorTheme: { name: 'uniform', params: { value: 3 } },
    };
    const freshJson = JSON.stringify(fresh);
    const shared = shareEqualParts(fresh, old);
    expect(JSON.stringify(shared)).toBe(freshJson);
    expect(JSON.stringify(fresh)).toBe(freshJson);
    expect(shared.colorTheme).toBe(old.colorTheme);
    expect(shared.type.params.x).toBe(old.type.params.x);
    expect(shared.type.params).not.toBe(old.type.params);
    expect(shareEqualParts({ ...old }, old)).toBe(old);
  });
});
