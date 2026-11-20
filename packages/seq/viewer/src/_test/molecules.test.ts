/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Molecules in the viewer (_spec/README.md "HoverModel", Molecules; "Behaviors", Links hover; "Commands", block.form):
// the hover's words name a modified residue's component and its alternatives, a molecule's residue, an atom and its
// residue; the hover carries the molecule's identifiers; `links.hover` lights a glycan's site from the glycan and the
// glycan from its site; `block.form` switches a block's form in one patch, a derived block written out first.

import type { Node } from '@molstar/core/repr/tree';
import { createSeqSpec, identifiers, modified } from '@molstar/seqx-core/lang';
import { afterAll, describe, expect, it } from 'vitest';
import { LinksHoverChannel, LinksHoverProvider } from '../behaviors/links-hover.js';
import { DefaultSeqViewerSpec } from '../default-spec.js';
import { Describer } from '../describe.js';
import { createSeqViewer, type SeqViewer } from '../viewer.js';

/** An Fc chain twice, its glycan on N10 of A, ATP with atoms; a section of the complex; a modified tRNA. */
function spec(): Node {
  const b = createSeqSpec({ title: 'molecules' });
  const fc = b.data.sequence('fc', {
    alphabet: 'protein',
    letters: 'MKTAYIAKQNQISFVKSHFSRQLEERLG',
    modified: modified([{ position: 0, component: 'MSE' }]),
  });
  const glycan = b.data.molecule('glycan', {
    class: 'branched',
    name: 'Fc glycan',
    residues: ['NAG', 'NAG', 'BMA'],
    identifiers: identifiers({ glytoucan: 'G15407YE' }),
  });
  const atp = b.data.molecule('atp', {
    residues: ['ATP'],
    atoms: { residue: [0, 0, 0], name: ['PG', 'O1G', 'N1'], element: ['P', 'O', 'N'] },
    identifiers: identifiers([['ccd', 'ATP']]),
  });
  const trna = b.data.sequence('trna', {
    alphabet: 'rna',
    letters: 'GCGGAUUUAGCUCAG',
    modified: modified([
      { position: 9, component: 'PSU' },
      { position: 12, component: 'H2U' },
      { position: 12, component: 'U' },
    ]),
  });
  b.data.compound('PSU', { name: 'pseudouridine' });
  b.data.complex('c', {}, (c) => {
    c.unit('A', { of: fc });
    c.unit('B', { of: fc });
    c.unit('g', { of: glycan });
    c.unit('l', { of: atp });
    c.links('links', {
      rows: { a_unit: ['A'], a_start: [9], a_end: [10], b_unit: ['g'], b_start: [0], b_end: [1], type: ['covalent'] },
    });
  });
  b.views.view('v', (v) => {
    v.section('complex', { of: '/data/c' }, (s) => s.track('letters').layer('letters', { type: 'letters' }));
    v.section('trna', { of: trna }, (s) => s.track('t').layer('letters', { type: 'letters' }));
  });
  return b.build();
}

const viewers: SeqViewer[] = [];
afterAll(() => {
  for (const v of viewers.splice(0)) v.retire();
});

async function make(links = true) {
  const base = DefaultSeqViewerSpec();
  const v = createSeqViewer({
    ...base,
    ui: undefined,
    behaviors: links ? [...base.behaviors, { provider: LinksHoverProvider, params: undefined }] : base.behaviors,
  });
  viewers.push(v);
  await v.load(spec());
  v.resize(800);
  await v.settle();
  return v;
}

const describeAt = (v: SeqViewer, section: string, system: string, position: number) => {
  const owner = {};
  v.hover.set({ section, system, position }, owner);
  const card = new Describer({
    spec: v.spec,
    store: v.store,
    addresses: v.addresses,
    view: v.view,
    layout: () => v.frame?.layout,
  }).describe(v.hover.current);
  const molecule = v.hover.molecule();
  v.hover.clear(owner);
  return {
    rows: Object.fromEntries((card?.rows ?? []).map((r) => [r.label, r.value])),
    text: card?.text ?? '',
    molecule,
  };
};

describe('molecules in the viewer', async () => {
  const v = await make();

  it('the hover names a modified residue’s component and name, and alternatives', () => {
    const psu = describeAt(v, '/views/v/trna', '/data/trna', 9);
    expect(psu.rows.Component).toBe('PSU (pseudouridine)');
    expect(psu.text).toContain('component PSU (pseudouridine)');
    expect(describeAt(v, '/views/v/trna', '/data/trna', 12).rows.Components).toBe('H2U | U');
    // a standard position says none
    expect(describeAt(v, '/views/v/trna', '/data/trna', 0).rows.Component).toBeUndefined();
    expect(describeAt(v, '/views/v/complex', '/data/c/A', 0).rows.Component).toBe('MSE');
  });

  it('a molecule’s residue by its code, an atom by its name and element with its residue; the molecule’s identifiers', () => {
    const nag = describeAt(v, '/views/v/complex', '/data/c/g', 2);
    expect(nag.rows.Component).toBe('BMA');
    expect(nag.molecule).toMatchObject({
      address: '/data/glycan',
      name: 'Fc glycan',
      class: 'branched',
      unit: '/data/c/g',
      identifiers: [{ namespace: 'glytoucan', value: 'G15407YE' }],
    });
    const atom = describeAt(v, '/views/v/complex', '/data/c/l/(atoms)', 1);
    expect([atom.rows.Atom, atom.rows.Residue]).toEqual(['O1G (O)', 'ATP 1']);
    expect(atom.molecule?.identifiers).toEqual([{ namespace: 'ccd', value: 'ATP' }]);
    expect(describeAt(v, '/views/v/trna', '/data/trna', 0).molecule).toBeUndefined();
  });

  it('links.hover: the glycan lights its site; the site lights the glycan; a hover elsewhere none', () => {
    const owner = {};
    v.hover.set({ section: '/views/v/complex', system: '/data/c/g', position: 1 }, owner);
    const site = v.highlight.get(LinksHoverChannel)?.location.parts.map((p) => [p.system, Array.from(p.ranges)]);
    expect(site).toEqual([['/data/c/A', [9, 10]]]);
    v.hover.set({ section: '/views/v/complex', system: '/data/c/A', position: 9 }, owner);
    const glycan = v.highlight.get(LinksHoverChannel)?.location.parts.map((p) => [p.system, Array.from(p.ranges)]);
    expect(glycan).toEqual([['/data/c/g', [0, 3]]]);
    v.hover.set({ section: '/views/v/complex', system: '/data/c/A', position: 3 }, owner);
    expect(v.highlight.get(LinksHoverChannel)).toBeUndefined();
    v.hover.clear(owner);
  });

  it('block.form: a derived block written out with its form, one undo step; the view otherwise the same', async () => {
    const w = await make(false);
    const before = w.view.resolved.value!.sections[0].blocks!.map((b) => [b.label, b.form]);
    expect(before).toEqual([
      ['A ×2', 'residues'],
      ['Fc glycan', 'whole'],
      ['ATP', 'whole'],
    ]);
    await w.execute('block.form', { address: '/views/v/complex/g', form: 'residues' });
    await w.settle();
    const after = w.view.resolved.value!.sections[0].blocks!.map((b) => [b.label, b.form]);
    expect(after).toEqual([
      ['A ×2', 'residues'],
      ['Fc glycan', 'residues'],
      ['ATP', 'whole'],
    ]);
    await w.execute('block.form', { address: '/views/v/complex/l', form: 'atoms' });
    expect(w.view.resolved.value!.sections[0].blocks![2].axis.system).toBe('/data/c/l/(atoms)');
  });
});
