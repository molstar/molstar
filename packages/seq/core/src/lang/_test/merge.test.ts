/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Merges over SeqViewSpec (`repr/patch`'s `merge3`): two 10-op edits of the p53 spec merged, the result consumed by
// `applyPatch` and by `diff`, two follows making a cycle refused, and the performance representative.

import {
  applyPatch,
  diff,
  type MergeOptions,
  merge3,
  type Op,
  type PatchValidator,
  patchFromNode,
  patchToNode,
} from '@molstar/core/repr/patch';
import { hashNode, type Node, resolveAddress } from '@molstar/core/repr/tree';
import { itWithin } from '@molstar/test-harness/timing';
import { describe, expect, it } from 'vitest';
import { checkSeqLang } from '../check.js';
import { SeqLanguage } from '../language.js';
import { p53Spec } from './p53.js';

const options: MergeOptions = { language: SeqLanguage };

/** Ten ops on the p53 spec: sets, an unset, an insert of a track with its layer, a reorder, a rename, a new sequence. */
const Ours: Op[] = [
  { op: 'set', at: [], key: 'title', value: 'p53, edited' },
  { op: 'set', at: ['data', 'tp53'], key: 'name', value: 'TP53' },
  { op: 'set', at: ['views', 'main', 'protein', 'letters'], key: 'height', value: 24 },
  {
    op: 'insert',
    parent: ['views', 'main', 'protein'],
    place: { last: true },
    node: {
      kind: 'track',
      id: 'conservation',
      children: [{ kind: 'layer', id: 'values', params: { type: 'values', data: '/data/conservation' } }],
    },
  },
  { op: 'move', at: ['views', 'main', 'gene'], parent: ['views', 'main'], place: { after: null } },
  { op: 'set', at: ['data', 'p53'], key: 'cyclic', value: false },
  {
    op: 'move',
    at: ['views', 'main', 'protein', 'letters'],
    parent: ['views', 'main', 'protein'],
    place: { after: null },
    rename: 'seq',
  },
  { op: 'set', at: ['views', 'main', 'protein'], key: 'range', value: [0, 20] },
  { op: 'unset', at: ['views', 'main', 'protein', 'domains'], key: 'label' },
  {
    op: 'insert',
    parent: ['data'],
    place: { last: true },
    node: { kind: 'sequence', id: 'p63', params: { alphabet: 'protein', letters: 'MNF' } },
  },
];

const track = (id: string, type: string, data?: string): Node => ({
  kind: 'track',
  id,
  children: [{ kind: 'layer', id: type, params: { type, ...(data ? { data } : {}) } }],
});

/** Ten other ops on it, none on what `Ours` changes: one insert placed after another. */
const Theirs: Op[] = [
  { op: 'set', at: ['views', 'main', 'protein', 'domains', 'ranges'], key: 'hidden', value: true },
  { op: 'set', at: ['views', 'main', 'protein', 'domains'], key: 'height', value: 30 },
  {
    op: 'insert',
    parent: ['views', 'main', 'gene'],
    place: { last: true },
    node: track('notes', 'values', '/data/conservation'),
  },
  { op: 'set', at: ['views', 'main', 'gene'], key: 'range', value: [5, 25] },
  {
    op: 'insert',
    parent: ['data'],
    place: { after: null },
    node: { kind: 'sequence', id: 'p73', params: { alphabet: 'protein', letters: 'MAQ' } },
  },
  { op: 'set', at: ['views', 'main', 'protein', 'domains', 'ranges'], key: 'color', value: 'red' },
  {
    op: 'insert',
    parent: ['views', 'main', 'protein'],
    place: { after: 'domains' },
    node: track('extra', 'letters'),
  },
  { op: 'set', at: ['data', 'p53'], key: 'name', value: 'P53' },
  { op: 'set', at: ['views', 'main'], key: 'layout', value: 'tracks' },
  { op: 'insert', parent: ['views', 'main', 'gene'], place: { after: 'notes' }, node: track('more', 'letters') },
];

function edited(spec: Node, ops: readonly Op[]): Node {
  const r = applyPatch(spec, { ops }, options);
  if (!r.ok) throw new Error(r.issues.map((i) => i.code).join());
  return r.root;
}

describe('merges over SeqViewSpec', () => {
  const base = p53Spec();
  const ours = edited(base, Ours);
  const theirs = edited(base, Theirs);

  it('merge two 10-op edits of the p53 spec as both applied in turn', () => {
    const r = merge3(base, ours, theirs, options);
    expect(r.ok).toBe(true);
    if (!r.ok) return;
    expect(r.issues).toEqual([]);
    // ours' track and theirs' both after `domains`: ours' first (siblings one side inserted, in either order)
    const inTurn = edited(ours, [
      ...Theirs.slice(0, 6),
      { ...(Theirs[6] as Op & { op: 'insert' }), place: { after: 'conservation' } },
      ...Theirs.slice(7),
    ]);
    expect(hashNode(r.root)).toBe(hashNode(inTurn));
    const ids = (n: Node | undefined) => n?.children?.map((c) => c.id);
    expect(ids(resolveAddress(r.root, ['views', 'main', 'protein']))).toEqual([
      'seq',
      'domains',
      'conservation',
      'extra',
    ]);
    const s = merge3(base, theirs, ours, options);
    expect(s.ok && ids(resolveAddress(s.root, ['views', 'main', 'protein']))).toEqual([
      'seq',
      'domains',
      'extra',
      'conservation',
    ]);
  });

  it('give a result consumed by applyPatch (its changes and root) and by diff (the rebased patch read back)', () => {
    const r = merge3(base, ours, theirs, options);
    if (!r.ok) throw new Error('conflicts');
    const forward = diff(ours, r.root);
    const again = applyPatch(ours, forward, { ...options, maxOps: Number.POSITIVE_INFINITY });
    expect(again.ok && hashNode(again.root)).toBe(hashNode(r.root));
    expect(again.ok && again.changes).toEqual(r.changes);
    const read = patchFromNode(patchToNode(forward));
    expect(read.issues).toEqual([]);
    const back = applyPatch(ours, read.patch!, { ...options, maxOps: Number.POSITIVE_INFINITY });
    expect(back.ok && hashNode(back.root)).toBe(hashNode(r.root));
  });

  it('refuse two follows that make a cycle, the partial ours', () => {
    const follows: PatchValidator = (_before, after) =>
      checkSeqLang(after).filter((i) => i.code === 'seq.lang.section.follow');
    const gene = ['views', 'main', 'gene'];
    const open = edited(base, [{ op: 'unset', at: gene, key: 'follow' }]);
    const x = edited(open, [{ op: 'set', at: gene, key: 'follow', value: '/views/main/protein' }]);
    const y = edited(open, [{ op: 'set', at: ['views', 'main', 'protein'], key: 'follow', value: '/views/main/gene' }]);
    const r = merge3(open, x, y, { ...options, validators: [follows] });
    expect(r.ok).toBe(false);
    if (r.ok) return;
    expect(r.conflicts.map((c) => c.kind)).toEqual(['invalid']);
    expect(resolveAddress(open, r.conflicts[0].at)?.kind).toBe('section');
    expect(r.conflicts[0].theirs).toEqual(diff(open, y).ops);
    expect(r.partial).toBe(x);
  });
});

describe('merges over SeqViewSpec: the performance representative', () => {
  const base = p53Spec();
  const ours = edited(base, Ours);
  const theirs = edited(base, Theirs);

  itWithin('a merge of two 10-op edits of the p53 spec', (within) => {
    merge3(base, ours, theirs, options); // warm
    // About 0.2 ms (two diffs and one apply); one sample, so 3 ms leaves room for a collection, under the 16 ms edit.
    const r = within(3, () => merge3(base, ours, theirs, options));
    expect(r.ok).toBe(true);
  });
});
