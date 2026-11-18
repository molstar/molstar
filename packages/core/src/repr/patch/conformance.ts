/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The patches' conformance suite (_spec/README.md "Patches", Tests): cases as data, each with its options, and
// `checkPatches`, which says how an `apply` breaks them. No test framework, as `repr/codec`'s suite: a language's own
// test runs it with its own cases (a Seq spec's), and expects no failure. Every op and its inverse at the root, in the
// middle and at the leaves; every refusal by code; normalization, refs following a move, rebinding and dangling refs
// refused and allowed, a remove-then-insert passing link integrity, validation of what changed, and a stored patch.

import { defineLanguage } from '../schema/check.js';
import type { Language } from '../schema/language.js';
import { hashNode, type Node, type Value } from '../tree/index.js';
import { patchLanguage, patchToNode } from './language.js';
import type { ApplyOptions, Changes, Patch, PatchResult } from './patch.js';

export interface PatchCase {
  readonly name: string;
  readonly root: Node;
  readonly patch: Patch;
  readonly options?: ApplyOptions;
  /** `codes`: the exact set of the result's issue codes. */
  readonly expect: {
    readonly ok: boolean;
    readonly root?: Node;
    /** The exact set of the result's issue codes. */
    readonly codes?: readonly string[];
    /** The result's `changes`, exactly. */
    readonly changes?: Changes;
    /** The result's `inverseNeeds`, exactly (none: absent). */
    readonly inverseNeeds?: { readonly allowDangling?: true; readonly allowRebind?: true };
  };
}

/** The language the cases are written in: boxes of leaves, a ref, a list of refs, a column, and a stored patch. */
export const PatchConformanceLanguage: Language = defineLanguage({
  id: 'pc',
  version: '1.0.0',
  root: 'doc',
  doc: "The patches' conformance trees.",
  kinds: {
    doc: { doc: 'The root.', parents: [], params: { title: { type: 'string', doc: 'A title.' } } },
    box: {
      doc: 'A box.',
      parents: ['doc', 'box'],
      ids: 'required',
      params: {
        label: { type: 'string', doc: 'A label.' },
        size: { type: 'number', min: 0, doc: 'A size.' },
        cols: { type: 'column', of: 'float32', doc: 'A column.' },
        link: { type: 'ref', doc: 'A ref.' },
        links: { type: 'list', of: { type: 'ref', doc: 'A ref.' }, doc: 'Refs.' },
        step: { type: 'node', language: 'pc.patch', doc: 'A stored patch.' },
        extra: { type: 'any', doc: 'Anything: a ref in it is text.' },
      },
    },
    leaf: {
      doc: 'A leaf.',
      parents: ['box'],
      ids: 'required',
      params: { label: { type: 'string', doc: 'A label.' }, link: { type: 'ref', doc: 'A ref.' } },
    },
    crate: { doc: 'A box of no leaves.', parents: ['doc'], ids: 'required', params: {} },
    bag: { doc: 'A bag of items without ids.', parents: ['doc'] },
    item: { doc: 'An item.', parents: ['bag'] },
  },
});

const languages = (id: string) => (id === 'pc.patch' ? patchLanguage(PatchConformanceLanguage) : undefined);
const lang: ApplyOptions = { language: PatchConformanceLanguage, validate: { languages } };

type P = Record<string, Value>;
const leaf = (id: string, params?: P): Node => ({ kind: 'leaf', id, ...(params ? { params } : {}) });
const box = (id: string, params: P | undefined, ...children: Node[]): Node => ({
  kind: 'box',
  id,
  ...(params ? { params } : {}),
  ...(children.length ? { children } : {}),
});
const doc = (params: P | undefined, ...children: Node[]): Node => ({
  kind: 'doc',
  ...(params ? { params } : {}),
  children,
});

const A = (params: P = { label: 'A', link: '/b/x' }, ...kids: Node[]) =>
  box('a', params, ...(kids.length ? kids : [leaf('x', { label: 'ax' }), leaf('y')]));
const B = (params: P = { label: 'B', link: '/a/x' }, ...kids: Node[]) =>
  box('b', params, ...(kids.length ? kids : [leaf('x')]));
const base = () => doc({ title: 'T' }, A(), B());
const bag: Node = { kind: 'bag', id: 'g', children: [{ kind: 'item' }, { kind: 'item' }] };
const deep = (n: number): Node => (n === 0 ? leaf('z') : box('z', undefined, deep(n - 1)));

const ok = (name: string, patch: Patch, root: Node, options?: ApplyOptions, codes?: readonly string[]): PatchCase => ({
  name,
  root: base(),
  patch,
  ...(options ? { options } : {}),
  expect: { ok: true, root, ...(codes ? { codes } : {}) },
});
const no = (name: string, patch: Patch, codes: readonly string[], options?: ApplyOptions, root?: Node): PatchCase => ({
  name,
  root: root ?? base(),
  patch,
  ...(options ? { options } : {}),
  expect: { ok: false, codes },
});

const step = (ops: Patch['ops']) => patchToNode({ ops });

/** The cases, each built fresh. */
function patchCases(): PatchCase[] {
  const T = 'f'.repeat(32);
  return [
    // every op, at the root, in the middle and at the leaves
    ok('set at the root', { ops: [{ op: 'set', at: [], key: 'title', value: 'U' }] }, doc({ title: 'U' }, A(), B())),
    ok(
      'set in the middle',
      { ops: [{ op: 'set', at: ['a'], key: 'label', value: 'A2' }] },
      doc({ title: 'T' }, A({ label: 'A2', link: '/b/x' }), B()),
    ),
    ok(
      'set at a leaf, in custom data',
      { ops: [{ op: 'set', at: ['a', 'y'], key: 'note', value: [1, { k: 'v' }], where: 'custom' }] },
      doc(
        { title: 'T' },
        A(undefined, leaf('x', { label: 'ax' }), { ...leaf('y'), custom: { note: [1, { k: 'v' }] } }),
        B(),
      ),
    ),
    ok('unset at the root', { ops: [{ op: 'unset', at: [], key: 'title' }] }, doc(undefined, A(), B())),
    ok(
      'unset at a leaf',
      { ops: [{ op: 'unset', at: ['a', 'x'], key: 'label' }] },
      doc({ title: 'T' }, A(undefined, leaf('x'), leaf('y')), B()),
    ),
    ok('unset of an absent key: a no-op', { ops: [{ op: 'unset', at: ['b'], key: 'size' }] }, base()),
    ok(
      'insert first, after a sibling and last',
      {
        ops: [
          { op: 'insert', parent: ['a'], place: { after: null }, node: leaf('w') },
          { op: 'insert', parent: ['a'], place: { after: 'x' }, node: leaf('xx') },
          { op: 'insert', parent: [], place: { last: true }, node: box('c', undefined) },
        ],
      },
      doc(
        { title: 'T' },
        A(undefined, leaf('w'), leaf('x', { label: 'ax' }), leaf('xx'), leaf('y')),
        B(),
        box('c', undefined),
      ),
    ),
    ok(
      'insert under a leaf',
      { ops: [{ op: 'insert', parent: ['b', 'x'], place: { last: true }, node: leaf('deeper') }] },
      doc({ title: 'T' }, A(), B(undefined, { ...leaf('x'), children: [leaf('deeper')] })),
    ),
    ok('remove in the middle', { ops: [{ op: 'remove', at: ['a'] }] }, doc({ title: 'T' }, B())),
    ok(
      'remove a leaf',
      { ops: [{ op: 'remove', at: ['a', 'y'] }] },
      doc({ title: 'T' }, A(undefined, leaf('x', { label: 'ax' })), B()),
    ),
    ok(
      'move a leaf to another parent, first',
      { ops: [{ op: 'move', at: ['a', 'y'], parent: ['b'], place: { after: null } }] },
      doc({ title: 'T' }, A(undefined, leaf('x', { label: 'ax' })), B(undefined, leaf('y'), leaf('x'))),
    ),
    ok(
      'move with a rename',
      { ops: [{ op: 'move', at: ['a', 'x'], parent: ['b'], place: { last: true }, rename: 'x2' }] },
      doc({ title: 'T' }, A(undefined, leaf('y')), B(undefined, leaf('x'), leaf('x2', { label: 'ax' }))),
    ),
    ok(
      'a reorder within a parent',
      { ops: [{ op: 'move', at: ['b'], parent: [], place: { after: null } }] },
      doc({ title: 'T' }, B(), A()),
    ),
    ok('move to where it is: a no-op', { ops: [{ op: 'move', at: ['b'], parent: [], place: { after: 'a' } }] }, base()),
    ok('replace the root', { ops: [{ op: 'replace', at: [], node: doc(undefined) }] }, doc(undefined)),
    ok(
      'replace in the middle, its kind changed',
      { ops: [{ op: 'replace', at: ['b'], node: { kind: 'crate', id: 'b' } }] },
      doc({ title: 'T' }, A(), { kind: 'crate', id: 'b' }),
    ),
    ok(
      'replace a leaf',
      { ops: [{ op: 'replace', at: ['a', 'y'], node: leaf('y', { label: 'new' }) }] },
      doc({ title: 'T' }, A(undefined, leaf('x', { label: 'ax' }), leaf('y', { label: 'new' })), B()),
    ),
    ok('replace with an equal node: a no-op', { ops: [{ op: 'replace', at: ['b'], node: B() }] }, base()),
    {
      name: 'a patch already applied (its result is the root)',
      root: base(),
      patch: { ops: [{ op: 'remove', at: ['zz'] }], base: T, result: hashNode(base()) },
      expect: { ok: true, root: base() },
    },
    // every refusal, by code
    no('a missing target', { ops: [{ op: 'set', at: ['zz'], key: 'label', value: 1 }] }, ['repr.patch.target']),
    no('a missing parent', { ops: [{ op: 'insert', parent: ['zz'], place: { last: true }, node: leaf('q') }] }, [
      'repr.patch.target',
    ]),
    no('an id a sibling has', { ops: [{ op: 'insert', parent: ['a'], place: { last: true }, node: leaf('x') }] }, [
      'repr.patch.id',
    ]),
    no(
      'an inserted node without an id',
      { ops: [{ op: 'insert', parent: ['a'], place: { last: true }, node: { kind: 'leaf' } }] },
      ['repr.patch.id'],
    ),
    no(
      'a malformed id',
      { ops: [{ op: 'move', at: ['a', 'x'], parent: ['a'], place: { last: true }, rename: 'a b' }] },
      ['repr.patch.id'],
    ),
    no(
      'a rename to a taken id',
      { ops: [{ op: 'move', at: ['a', 'x'], parent: ['a'], place: { last: true }, rename: 'y' }] },
      ['repr.patch.id'],
    ),
    no('a replacement with another id', { ops: [{ op: 'replace', at: ['a', 'x'], node: leaf('q') }] }, [
      'repr.patch.id',
    ]),
    no(
      'a place naming no sibling',
      { ops: [{ op: 'insert', parent: ['a'], place: { after: 'zz' }, node: leaf('q') }] },
      ['repr.patch.place'],
    ),
    no(
      'a place after the moved node itself',
      { ops: [{ op: 'move', at: ['a', 'x'], parent: ['a'], place: { after: 'x' } }] },
      ['repr.patch.place'],
    ),
    no('a move under itself', { ops: [{ op: 'move', at: ['a'], parent: ['a', 'x'], place: { last: true } }] }, [
      'repr.patch.cycle',
    ]),
    no('a reserved key', { ops: [{ op: 'set', at: ['a'], key: 'kind', value: 'leaf' }] }, ['repr.patch.key']),
    no('a reserved $ key', { ops: [{ op: 'unset', at: ['a'], key: '$column' }] }, ['repr.patch.key']),
    no(
      'an insert under a parent with children without ids',
      { ops: [{ op: 'insert', parent: ['g'], place: { last: true }, node: { kind: 'item', id: 'i' } }] },
      ['repr.patch.unaddressable'],
      undefined,
      doc(undefined, bag),
    ),
    no('a remove of the root', { ops: [{ op: 'remove', at: [] }] }, ['repr.patch.op']),
    no('a move of the root', { ops: [{ op: 'move', at: [], parent: ['a'], place: { last: true } }] }, [
      'repr.patch.op',
    ]),
    no('an unknown op', { ops: [{ op: 'swap', at: ['a'] } as never] }, ['repr.patch.op']),
    no(
      'more ops than allowed',
      {
        ops: [
          { op: 'set', at: [], key: 'title', value: 'U' },
          { op: 'set', at: [], key: 'title', value: 'V' },
        ],
      },
      ['repr.patch.limit'],
      { maxOps: 1 },
    ),
    no(
      'a node past the depth bound',
      { ops: [{ op: 'insert', parent: ['a'], place: { last: true }, node: deep(600) }] },
      ['repr.patch.depth'],
    ),
    no('a stale patch', { ops: [], base: T }, ['repr.patch.stale']),
    no('a patch without a base, where one is required', { ops: [] }, ['repr.patch.base'], { requireBase: true }),
    no('a patch that gives another result', { ops: [{ op: 'set', at: [], key: 'title', value: 'U' }], result: T }, [
      'repr.patch.result',
    ]),
    no(
      'a later op refused refuses the whole patch',
      {
        ops: [
          { op: 'set', at: [], key: 'title', value: 'U' },
          { op: 'remove', at: ['zz'] },
        ],
      },
      ['repr.patch.target'],
    ),
    // with a language: normalization, refs, link integrity, validation
    ok(
      'a plain array set into a column is normalized',
      { ops: [{ op: 'set', at: ['b'], key: 'cols', value: [1, 2] }] },
      doc({ title: 'T' }, A(), B({ label: 'B', link: '/a/x', cols: new Float32Array([1, 2]) })),
      lang,
    ),
    ok(
      'a ref follows a rename, and back by the inverse',
      { ops: [{ op: 'move', at: ['a'], parent: [], place: { after: null }, rename: 'c' }] },
      doc(
        { title: 'T' },
        box('c', { label: 'A', link: '/b/x' }, leaf('x', { label: 'ax' }), leaf('y')),
        B({ label: 'B', link: '/c/x' }),
      ),
      lang,
    ),
    {
      name: 'refs in a list follow a move to another parent, a form after the id run kept',
      root: doc({ title: 'T' }, A({ label: 'A', link: '/b/x', links: ['/a/x/(range :from 1 :to 2)', '/a/y'] }), B()),
      patch: { ops: [{ op: 'move', at: ['a', 'x'], parent: ['b'], place: { last: true }, rename: 'ax' }] },
      options: lang,
      expect: {
        ok: true,
        root: doc(
          { title: 'T' },
          A({ label: 'A', link: '/b/x', links: ['/b/ax/(range :from 1 :to 2)', '/a/y'] }, leaf('y')),
          B({ label: 'B', link: '/b/ax' }, leaf('x'), leaf('ax', { label: 'ax' })),
        ),
      },
    },
    no('a dangling ref refused', { ops: [{ op: 'remove', at: ['b'] }] }, ['repr.patch.dangling'], lang),
    ok(
      'a dangling ref allowed, its inverse needing a rebinding',
      { ops: [{ op: 'remove', at: ['b'] }] },
      doc({ title: 'T' }, A()),
      { ...lang, allowDangling: true },
      ['repr.patch.dangling'],
    ),
    no(
      'a ref written dangling refused',
      { ops: [{ op: 'set', at: ['a'], key: 'link', value: '/nowhere' }] },
      ['repr.patch.dangling'],
      lang,
    ),
    {
      name: 'a rebinding refused',
      root: doc(undefined, box('a', { link: '/z/x' })),
      patch: { ops: [{ op: 'insert', parent: [], place: { last: true }, node: box('z', undefined, leaf('x')) }] },
      options: lang,
      expect: { ok: false, codes: ['repr.patch.rebound'] },
    },
    {
      name: 'a rebinding allowed, by a move: its inverse restores the ref the move back would take along',
      root: doc(undefined, box('a', { link: '/z/x' }), box('m', { link: '/m/x' }, leaf('x'))),
      patch: { ops: [{ op: 'move', at: ['m'], parent: [], place: { last: true }, rename: 'z' }] },
      options: { ...lang, allowRebind: true },
      expect: { ok: true, root: doc(undefined, box('a', { link: '/z/x' }), box('z', { link: '/z/x' }, leaf('x'))) },
    },
    ok(
      'a remove then an insert (a diff) below what a ref names passes link integrity',
      {
        ops: [
          { op: 'remove', at: ['a', 'y'] },
          { op: 'insert', parent: ['a'], place: { last: true }, node: leaf('y', { label: 'ay' }) },
          { op: 'set', at: ['b'], key: 'link', value: '/a' },
        ],
      },
      doc(
        { title: 'T' },
        A(undefined, leaf('x', { label: 'ax' }), leaf('y', { label: 'ay' })),
        B({ label: 'B', link: '/a' }),
      ),
      lang,
    ),
    no(
      'a remove then an insert at the address a ref names rebinds it: refused',
      {
        ops: [
          { op: 'remove', at: ['a', 'x'] },
          { op: 'insert', parent: ['a'], place: { after: null }, node: leaf('x', { label: 'ax2' }) },
        ],
      },
      ['repr.patch.rebound'],
      lang,
    ),
    ok(
      'a remove then an insert at the address a ref names, the rebinding allowed',
      {
        ops: [
          { op: 'remove', at: ['a', 'x'] },
          { op: 'insert', parent: ['a'], place: { after: null }, node: leaf('x', { label: 'ax2' }) },
        ],
      },
      doc({ title: 'T' }, A(undefined, leaf('x', { label: 'ax2' }), leaf('y')), B()),
      { ...lang, allowRebind: true },
    ),
    {
      name: 'an inserted and a replaced subtree are normalized',
      root: base(),
      patch: {
        ops: [
          {
            op: 'insert',
            parent: [],
            place: { last: true },
            node: box('c', { cols: [1, 2] }, box('d', { cols: [3] })),
          },
          { op: 'replace', at: ['b'], node: box('b', { label: 'B', link: '/a/x', cols: [4] }, leaf('x')) },
        ],
      },
      options: lang,
      expect: {
        ok: true,
        root: doc(
          { title: 'T' },
          A(),
          box('b', { label: 'B', link: '/a/x', cols: new Float32Array([4]) }, leaf('x')),
          box('c', { cols: new Float32Array([1, 2]) }, box('d', { cols: new Float32Array([3]) })),
        ),
        changes: { made: [['b'], ['c']], params: [], children: [[]], removed: [], moved: [] },
      },
    },
    {
      name: 'a ref follows by its id run, not by its text: /ab stays when /a moves',
      root: doc(
        undefined,
        box('a', undefined, leaf('x')),
        box('ab', { link: '/ab', links: ['/a/x', '/ab/(range :from 1)'] }),
      ),
      patch: { ops: [{ op: 'move', at: ['a'], parent: [], place: { last: true }, rename: 'c' }] },
      options: lang,
      expect: {
        ok: true,
        root: doc(
          undefined,
          box('ab', { link: '/ab', links: ['/c/x', '/ab/(range :from 1)'] }),
          box('c', undefined, leaf('x')),
        ),
        changes: {
          made: [],
          params: [{ at: ['ab'], keys: ['links'], where: 'params' }],
          children: [[]],
          removed: [],
          moved: [{ from: ['a'], to: ['c'] }],
        },
      },
    },
    {
      name: 'refs in custom data and in an any param are not rewritten',
      root: doc(undefined, box('a', undefined, leaf('x')), {
        ...box('r', { extra: '/a/x' }),
        custom: { link: '/a/x' },
      }),
      patch: { ops: [{ op: 'move', at: ['a'], parent: [], place: { last: true }, rename: 'c' }] },
      options: lang,
      expect: {
        ok: true,
        root: doc(
          undefined,
          { ...box('r', { extra: '/a/x' }), custom: { link: '/a/x' } },
          box('c', undefined, leaf('x')),
        ),
      },
    },
    {
      name: 'refs inside the moved subtree to itself move with it',
      root: doc(undefined, box('a', { link: '/a/x' }, leaf('x', { link: '/a' })), box('b', undefined)),
      patch: { ops: [{ op: 'move', at: ['a'], parent: ['b'], place: { last: true } }] },
      options: lang,
      expect: {
        ok: true,
        root: doc(undefined, box('b', undefined, box('a', { link: '/b/a/x' }, leaf('x', { link: '/b/a' })))),
        changes: {
          made: [],
          params: [
            { at: ['b', 'a'], keys: ['link'], where: 'params' },
            { at: ['b', 'a', 'x'], keys: ['link'], where: 'params' },
          ],
          children: [[], ['b']],
          removed: [],
          moved: [{ from: ['a'], to: ['b', 'a'] }],
        },
      },
    },
    {
      name: 'ops that interact: a ref written then moved along, an inserted node moved, a set at its new address',
      root: doc(undefined, box('a', undefined, leaf('x')), box('b', undefined)),
      patch: {
        ops: [
          { op: 'set', at: ['b'], key: 'link', value: '/a/x' },
          { op: 'insert', parent: ['a'], place: { last: true }, node: leaf('n', { link: '/a/x' }) },
          { op: 'move', at: ['a'], parent: ['b'], place: { last: true }, rename: 'c' },
          { op: 'set', at: ['b', 'c', 'n'], key: 'label', value: 'N' },
          { op: 'move', at: ['b', 'c', 'x'], parent: ['b', 'c'], place: { last: true }, rename: 'x2' },
        ],
      },
      options: lang,
      expect: {
        ok: true,
        root: doc(
          undefined,
          box('b', { link: '/b/c/x2' }, box('c', undefined, leaf('n', { link: '/b/c/x2', label: 'N' }), leaf('x2'))),
        ),
      },
    },
    {
      name: 'a dangling ref allowed: its inverse needs a rebinding, and says so',
      root: base(),
      patch: { ops: [{ op: 'remove', at: ['b'] }] },
      options: { ...lang, allowDangling: true },
      expect: {
        ok: true,
        root: doc({ title: 'T' }, A()),
        codes: ['repr.patch.dangling'],
        inverseNeeds: { allowRebind: true },
        changes: { made: [], params: [], children: [[]], removed: [['b']], moved: [] },
      },
    },
    {
      name: 'a remove of content holding a ref that already dangled: the inverse restores it',
      root: doc(undefined, box('a', { link: '/zz' }, leaf('x', { link: '/zz/y' }))),
      patch: { ops: [{ op: 'remove', at: ['a'] }] },
      options: lang,
      expect: { ok: true, root: doc(undefined) },
    },
    {
      name: 'the inverse of a move never drags a ref the move did not rewrite',
      root: doc(undefined, box('a', { link: '/c/x' }), box('m', undefined, leaf('x'))),
      patch: { ops: [{ op: 'move', at: ['m'], parent: [], place: { last: true }, rename: 'c' }] },
      options: { ...lang, allowRebind: true },
      expect: { ok: true, root: doc(undefined, box('a', { link: '/c/x' }), box('c', undefined, leaf('x'))) },
    },
    no(
      'a value outside its range refused',
      { ops: [{ op: 'set', at: ['a'], key: 'size', value: -1 }] },
      ['repr.schema.param.range'],
      lang,
    ),
    no(
      'a node under a parent it may not sit under refused',
      { ops: [{ op: 'insert', parent: [], place: { last: true }, node: leaf('q') }] },
      ['repr.schema.parent.invalid'],
      lang,
    ),
    no(
      'a kind-changing replace validates its kept children',
      { ops: [{ op: 'replace', at: ['a'], node: { kind: 'crate', id: 'a', children: [leaf('x'), leaf('y')] } }] },
      ['repr.schema.parent.invalid'],
      lang,
    ),
    {
      name: 'an error present before is not reported again by an edit below it',
      root: doc(undefined, box('a', { size: -5 }, leaf('x'))),
      patch: { ops: [{ op: 'set', at: ['a', 'x'], key: 'label', value: 'ok' }] },
      options: lang,
      expect: { ok: true, root: doc(undefined, box('a', { size: -5 }, leaf('x', { label: 'ok' }))) },
    },
    {
      name: 'a stored patch inserts a node and sets a key on it; both follow a rename of their parent',
      root: doc(
        undefined,
        box('a', undefined, leaf('x')),
        box('s', {
          step: step([
            { op: 'insert', parent: ['a'], place: { last: true }, node: leaf('n') },
            { op: 'set', at: ['a', 'n'], key: 'label', value: 'new' },
          ]),
        }),
      ),
      patch: { ops: [{ op: 'move', at: ['a'], parent: [], place: { after: null }, rename: 'c' }] },
      options: lang,
      expect: {
        ok: true,
        root: doc(
          undefined,
          box('c', undefined, leaf('x')),
          box('s', {
            step: step([
              { op: 'insert', parent: ['c'], place: { last: true }, node: leaf('n') },
              { op: 'set', at: ['c', 'n'], key: 'label', value: 'new' },
            ]),
          }),
        ),
      },
    },
    {
      name: 'a stored set of a ref value is rewritten by a rename of what it names',
      root: doc(
        undefined,
        box('a', undefined, leaf('x')),
        box('s', { step: step([{ op: 'set', at: ['s'], key: 'link', value: '/a/x' }]) }),
      ),
      patch: { ops: [{ op: 'move', at: ['a', 'x'], parent: ['a'], place: { last: true }, rename: 'x2' }] },
      options: lang,
      expect: {
        ok: true,
        root: doc(
          undefined,
          box('a', undefined, leaf('x2')),
          box('s', { step: step([{ op: 'set', at: ['s'], key: 'link', value: '/a/x2' }]) }),
        ),
      },
    },
    {
      name: 'a stored patch counts in link integrity',
      root: doc(
        undefined,
        box('a', undefined, leaf('x')),
        box('s', { step: step([{ op: 'set', at: ['a', 'x'], key: 'label', value: 'v' }]) }),
      ),
      patch: { ops: [{ op: 'remove', at: ['a', 'x'] }] },
      options: lang,
      expect: { ok: false, codes: ['repr.patch.dangling'] },
    },
  ];
}

/** The cases. */
export const PatchCases: readonly PatchCase[] = patchCases();

/**
 * How `apply` breaks the cases (empty when it holds them all): a result other than expected, a root other than
 * expected (by hash), issue codes other than expected, the caller's root changed, or an inverse that does not give
 * the root back when applied with what it needs.
 */
export function checkPatches(
  apply: (root: Node, patch: Patch, options?: ApplyOptions) => PatchResult,
  cases: readonly PatchCase[] = PatchCases,
): readonly string[] {
  const failures: string[] = [];
  for (const c of cases) {
    const fail = (why: string) => failures.push(`${c.name}: ${why}`);
    const before = hashNode(structuredClone(c.root));
    let r: PatchResult;
    try {
      r = apply(c.root, c.patch, c.options);
    } catch (e) {
      fail(`apply threw (${e instanceof Error ? e.message : String(e)})`);
      continue;
    }
    if (hashNode(structuredClone(c.root)) !== before) fail("the caller's root changed");
    if (r.ok !== c.expect.ok)
      fail(
        `expected ${c.expect.ok ? 'it applied' : 'a refusal'}, got ${r.ok ? 'it applied' : r.issues.map((i) => i.code).join(', ')}`,
      );
    if (c.expect.codes) {
      const got = [...new Set(r.issues.map((i) => i.code))].sort();
      const want = [...new Set(c.expect.codes)].sort();
      if (got.join() !== want.join())
        fail(`expected the codes ${want.join(', ') || 'none'}, got ${got.join(', ') || 'none'}`);
    }
    if (!r.ok) continue;
    if (c.expect.root && hashNode(r.root) !== hashNode(c.expect.root)) fail('the root is not the expected one');
    if (r.revision !== hashNode(r.root)) fail("the revision is not the root's hash");
    if (!r.changed && hashNode(r.root) !== hashNode(c.root)) fail('unchanged, and the root changed');
    if (r.changed && !r.inverse) fail('changed, and no inverse');
    if (!r.inverse) continue;
    if (c.expect.changes && JSON.stringify(r.changes) !== JSON.stringify(c.expect.changes))
      fail(`the changes are ${JSON.stringify(r.changes)}`);
    if ('inverseNeeds' in c.expect && JSON.stringify(r.inverseNeeds) !== JSON.stringify(c.expect.inverseNeeds))
      fail(`the inverse needs ${JSON.stringify(r.inverseNeeds ?? {})}`);
    // the inverse with what it needs only: the case's own allowances are not its
    const { allowDangling: _d, allowRebind: _r, ...own } = c.options ?? {};
    const back = apply(r.root, r.inverse, { ...own, ...r.inverseNeeds, maxOps: Number.POSITIVE_INFINITY });
    if (!back.ok) fail(`the inverse was refused (${back.issues.map((i) => i.code).join(', ')})`);
    else if (hashNode(back.root) !== hashNode(c.root)) fail('the inverse does not give the root back');
  }
  return failures;
}
