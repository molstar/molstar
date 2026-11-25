/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Patches over SeqViewSpec: `repr/patch`'s conformance suite run with Seq's own cases over `SeqLanguage` (every child
// has an id, so every parent takes inserts, removes and moves), and the performance representative, a 10-op patch on
// the p53 spec applied and diffed back. Then the cases with the checks of a change and the story steps' slots, as the
// viewer's spec model applies them.

import { applyPatch, diff, type Op, type Patch, patchToNode } from '@molstar/core/repr/patch';
import { checkPatches, type PatchCase, PatchCases } from '@molstar/core/repr/patch/conformance';
import { hashNode, type Node } from '@molstar/core/repr/tree';
import { itWithin } from '@molstar/test-harness/timing';
import { describe, expect, it } from 'vitest';
import { createSeqSpec } from '../builder.js';
import { SeqPatchValidator } from '../check-changes.js';
import { SeqLanguage, SeqValidateOptions } from '../language.js';
import { P53_LETTERS, p53Spec } from './p53.js';

const options = { language: SeqLanguage };

/** The p53 spec with its protein section renamed and its domains track under the gene section. */
function p53Moved(): Node {
  const b = createSeqSpec({ title: 'p53' });
  const p53 = b.data.sequence('p53', { alphabet: 'protein', letters: P53_LETTERS });
  const chr17 = b.data.source('chr17', {
    format: 'fasta',
    uri: 'data/chr17-region.fa',
    index: 'data/chr17-region.fa.fai',
  });
  const gene = b.data.sequence('tp53', { alphabet: 'dna', source: chr17 });
  b.data.map('p53-on-gene', {
    from: p53,
    to: gene,
    method: 'codon',
    exons: { start: [100, 160], end: [150, 200] },
    strand: '-',
  });
  const domains = b.data.ranges('domains', {
    target: p53,
    rows: { start: new Int32Array([0, 17]), end: new Int32Array([14, 30]), label: ['Disordered', 'TAD'] },
  });
  b.data.values('conservation', { target: p53, values: new Float32Array(P53_LETTERS.length).fill(0.5) });
  b.views.view('main', (v) => {
    const protein = v.section('prot', { of: p53 });
    protein.track('letters').layer('letters', { type: 'letters' });
    v.section('gene', { of: gene, follow: protein }, (g) =>
      g.track('domains', { label: 'Domains' }).layer('ranges', { type: 'ranges', data: domains }),
    );
  });
  return b.build();
}

/** `root` with the node at `address` changed by `f`. */
function edit(root: Node, address: readonly string[], f: (n: Node) => Node): Node {
  if (!address.length) return f(root);
  const [id, ...rest] = address;
  return { ...root, children: root.children!.map((c) => (c.id === id ? edit(c, rest, f) : c)) };
}

const SeqCases: PatchCase[] = [
  {
    name: 'a layer hidden: a one-key set on it',
    root: p53Spec(),
    patch: {
      ops: [
        {
          op: 'set',
          at: ['views', 'main', 'protein', 'letters', 'letters'],
          key: 'hidden',
          value: true,
        },
      ],
    },
    options,
    expect: {
      ok: true,
      root: edit(p53Spec(), ['views', 'main', 'protein', 'letters', 'letters'], (n) => ({
        ...n,
        params: { ...n.params, hidden: true },
      })),
    },
  },
  {
    name: 'a track moved to another section, and a section renamed with the follow ref following',
    root: p53Spec(),
    patch: {
      ops: [
        {
          op: 'move',
          at: ['views', 'main', 'protein', 'domains'],
          parent: ['views', 'main', 'gene'],
          place: { last: true },
        },
        {
          op: 'move',
          at: ['views', 'main', 'protein'],
          parent: ['views', 'main'],
          place: { after: null },
          rename: 'prot',
        },
      ],
    },
    options,
    expect: { ok: true, root: p53Moved() },
  },
  {
    name: 'a dangling ref refused: the data a layer draws removed',
    root: p53Spec(),
    patch: { ops: [{ op: 'remove', at: ['data', 'domains'] }] },
    options,
    expect: { ok: false, codes: ['repr.patch.dangling'] },
  },
  {
    name: 'a layer without an id refused by the language',
    root: p53Spec(),
    patch: {
      ops: [
        {
          op: 'replace',
          at: ['views', 'main', 'protein', 'letters'],
          node: { kind: 'track', id: 'letters', children: [{ kind: 'layer', params: { type: 'letters' } }] },
        },
      ],
    },
    options,
    expect: { ok: false, codes: ['repr.schema.id.missing'] },
  },
];

/** What the viewer's `SpecModel` applies with: the language, its slots, the checks of a change. */
const checked = { language: SeqLanguage, validate: SeqValidateOptions, validators: [SeqPatchValidator] };

/** p53 with a story step that sets the domains track's label and inserts a track after it. */
function p53WithStep(): Node {
  const step: Patch = {
    ops: [
      { op: 'set', at: ['views', 'main', 'protein', 'domains'], key: 'label', value: 'Domains (step)' },
      {
        op: 'insert',
        parent: ['views', 'main', 'protein'],
        place: { after: 'domains' },
        node: { kind: 'track', id: 'note', children: [{ kind: 'layer', id: 'n', params: { type: 'letters' } }] },
      },
    ],
  };
  return edit(p53Spec(), [], (n) => ({
    ...n,
    children: [
      ...n.children!,
      {
        kind: 'snapshots',
        id: 'snapshots',
        children: [{ kind: 'snapshot', id: 'one', params: { patch: patchToNode(step) } }],
      },
    ],
  }));
}

/** Seq's cases with the checks of a change and the stored patches' slots: what the viewer's spec model applies. */
const CheckedCases: PatchCase[] = [
  {
    name: 'a range past the protein refused by the checks of a change',
    root: p53Spec(),
    patch: { ops: [{ op: 'set', at: ['views', 'main', 'protein'], key: 'range', value: [0, 300] }] },
    options: checked,
    expect: { ok: false, codes: ['seq.lang.section.range'] },
  },
  {
    name: 'the protein shortened under its domains and values refused',
    root: p53Spec(),
    patch: { ops: [{ op: 'set', at: ['data', 'p53'], key: 'letters', value: 'MEEP' }] },
    options: checked,
    expect: { ok: false, codes: ['seq.lang.ranges.bounds', 'seq.lang.values.length'] },
  },
  {
    name: 'a step’s target renamed: the stored patch follows, and so does its inverse',
    root: p53WithStep(),
    patch: {
      ops: [
        {
          op: 'move',
          at: ['views', 'main', 'protein', 'domains'],
          parent: ['views', 'main', 'protein'],
          place: { last: true },
          rename: 'dom',
        },
      ],
    },
    options: checked,
    expect: { ok: true },
  },
  {
    name: 'a removal of what a step names refused: its addresses are refs of the spec',
    root: p53WithStep(),
    patch: { ops: [{ op: 'remove', at: ['views', 'main', 'protein', 'domains'] }] },
    options: checked,
    expect: { ok: false, codes: ['repr.patch.dangling'] },
  },
  {
    name: 'a step whose patch has a base refused by the slot',
    root: p53Spec(),
    patch: {
      ops: [
        {
          op: 'insert',
          parent: [],
          place: { last: true },
          node: {
            kind: 'snapshots',
            id: 'snapshots',
            children: [{ kind: 'snapshot', id: 's', params: { patch: { kind: 'patch', params: { base: 'x' } } } }],
          },
        },
      ],
    },
    options: checked,
    expect: { ok: false, codes: ['repr.schema.param.pattern'] },
  },
];

describe('patches over SeqViewSpec', () => {
  it("hold the core's cases and Seq's own", () => {
    expect(checkPatches(applyPatch, PatchCases)).toEqual([]);
    expect(checkPatches(applyPatch, SeqCases)).toEqual([]);
    expect(checkPatches(applyPatch, CheckedCases)).toEqual([]);
  });
});

/** Ten ops on the p53 spec: sets, an unset, an insert of a track with its layer, a reorder, a rename, a new sequence. */
const TenOps: Op[] = [
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

describe('patches over SeqViewSpec: the performance representative', () => {
  const spec = p53Spec();
  const roundTrip = () => {
    const forward = applyPatch(spec, { ops: TenOps }, options);
    if (!forward.ok) throw new Error(forward.issues.map((i) => i.code).join());
    const back: Patch = diff(forward.root, spec);
    const restored = applyPatch(forward.root, back, { ...options, maxOps: Number.POSITIVE_INFINITY });
    if (!restored.ok) throw new Error(restored.issues.map((i) => i.code).join());
    return { forward, back, restored };
  };

  itWithin('a 10-op patch on the p53 spec, applied and diffed back', (within) => {
    roundTrip(); // warm
    // About 0.11 ms; one sample, so 5 ms leaves room for a collection (1 ms read 1.2 under a parallel heavy run).
    const { forward, back, restored } = within(5, roundTrip);
    expect(forward.changed).toBe(true);
    expect(back.ops.length).toBeGreaterThan(0);
    expect(hashNode(restored.root)).toBe(hashNode(spec));
  });
});
