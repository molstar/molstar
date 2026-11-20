/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Molecules in the viewer (_spec/README.md "HoverModel", Molecules; "Behaviors", Links hover; "Commands", block.form):
// the hover's words name a modified residue's component and its alternatives, a molecule's residue, an atom and its
// residue; the hover carries the molecule's identifiers; `links.hover` lights a glycan's site from the glycan and the
// glycan from its site; `block.form` switches a block's form in one patch, a derived block by one `block` child that overrides it alone.

import type { Node } from '@molstar/core/repr/tree';
import { createSeqSpec, identifiers, modified } from '@molstar/seqx-core/lang';
import { afterAll, describe, expect, it } from 'vitest';
import { LinksHoverChannel, LinksHoverProvider } from '../behaviors/links-hover.js';
import { DefaultSeqViewerSpec } from '../default-spec.js';
import { Describer } from '../describe.js';
import { type UiHistory, uiHistory } from '../ui-history.js';
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
  // the links' sites as ranges: their ends plans built when this view resolves
  b.views.view('sites', (v) =>
    v.section('complex', { of: '/data/c' }, (s) =>
      s.track('sites').layer('sites', { type: 'ranges', data: '/data/c/links' as never }),
    ),
  );
  return b.build();
}

/** ATP twice, merged into one block by its atoms (an override of the first unit), a sodium bound to the second's O1G. */
function atomsSpec(): Node {
  const b = createSeqSpec({ title: 'atoms' });
  const atp = b.data.molecule('atp', {
    residues: ['ATP'],
    atoms: { residue: [0, 0, 0], name: ['PG', 'O1G', 'N1'], element: ['P', 'O', 'N'] },
  });
  const na = b.data.molecule('na', { class: 'ion', residues: ['NA'] });
  b.data.complex('c', {}, (c) => {
    c.unit('l', { of: atp });
    c.unit('l2', { of: atp });
    c.unit('n', { of: na });
    c.links('links', {
      rows: {
        a_unit: ['n'],
        a_start: [0],
        a_end: [1],
        b_unit: ['l2'],
        b_start: [0],
        b_end: [1],
        type: ['metal'],
        a_atom: ['NA'],
        b_atom: ['O1G'],
      },
    });
  });
  b.views.view('v', (v) =>
    v.section('complex', { of: '/data/c' }, (s) => {
      s.block('l', { form: 'atoms' });
      s.track('letters').layer('letters', { type: 'letters' });
    }),
  );
  return b.build();
}

const viewers: SeqViewer[] = [];
afterAll(() => {
  for (const v of viewers.splice(0)) v.retire();
});

async function make(links = true, history = false, root: Node = spec()) {
  const base = DefaultSeqViewerSpec();
  const v = createSeqViewer({
    ...base,
    ui: undefined,
    ...(history ? { history: uiHistory() } : {}),
    behaviors: links ? [...base.behaviors, { provider: LinksHoverProvider, params: undefined }] : base.behaviors,
  });
  viewers.push(v);
  await v.load(root);
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

  it('links.hover onto an atoms-form merged block: the copy’s residue lit as its atoms, on the block’s atoms system', async () => {
    const w = await make(true, false, atomsSpec());
    expect(w.view.resolved.value!.sections[0].blocks!.map((b) => [b.label, b.axis.system])).toEqual([
      ['ATP ×2', '/data/c/l/(atoms)'],
      ['NA', '/data/c/n'],
    ]);
    const owner = {};
    w.hover.set({ section: '/views/v/complex', system: '/data/c/n', position: 0 }, owner);
    const lit = w.highlight.get(LinksHoverChannel)?.location.parts.map((p) => [p.system, Array.from(p.ranges)]);
    expect(lit).toEqual([['/data/c/l/(atoms)', [0, 3]]]);
    // and back: an atom of the merged block lights the sodium
    w.hover.set({ section: '/views/v/complex', system: '/data/c/l/(atoms)', position: 1, block: 0 }, owner);
    const back = w.highlight.get(LinksHoverChannel)?.location.parts.map((p) => [p.system, Array.from(p.ranges)]);
    expect(back).toEqual([['/data/c/n', [0, 1]]]);
    w.hover.clear(owner);
  });

  it('a view shown later with a ranges layer over the links: its ends plans built then, and opened into the store', async () => {
    const w = await make(false);
    w.view.show('/views/v');
    // an edit while another view is shown: a model whose plans no view has asked for
    await w.spec.patch('Rename', { ops: [{ op: 'set', at: ['data', 'c'], key: 'name', value: 'Fc' }] });
    await w.settle();
    const ends = w.spec.model.value.annotations.get('/data/c/links')!.ends!;
    expect([...ends.built()]).toEqual([]);
    const opened = () => w.store.store.coverage('/data/c/links#/data/c/A') !== undefined;
    expect(opened()).toBe(false);
    w.view.show('/views/sites');
    await w.settle();
    expect(w.view.resolved.value!.sections[0].blocks![0].layers[0][0].placed).toBe(true);
    expect(opened()).toBe(true);
  });

  it('a click on a chip selects its every residue', async () => {
    v.draw();
    const p = v.frame!.layout.sections[0].placement;
    const i = p.strips.axes.findIndex((a) => a.system === '/data/c/g');
    const a = p.strips.axes[i];
    const hit = v.hit(a.x0 + a.width / 2, p.strips.bandY[i * p.strips.bands] + 2);
    expect(hit?.system).toBe('/data/c/g');
    await v.execute('selection.set', { hit });
    const parts = v.selection.value.value.shown.parts.filter((q) => q.system === '/data/c/g');
    expect(parts.map((q) => Array.from(q.ranges))).toEqual([[0, 3]]);
  });

  it('block.form: one `block` child overrides that unit, the others stay derived; one undo step; a unit added later is listed', async () => {
    const w = await make(false, true);
    const forms = () => w.view.resolved.value!.sections[0].blocks!.map((b) => [b.label, b.form]);
    expect(forms()).toEqual([
      ['A ×2', 'residues'],
      ['Fc glycan', 'whole'],
      ['ATP', 'whole'],
    ]);
    await w.execute('block.form', { address: '/views/v/complex/g', form: 'residues' });
    await w.settle();
    expect(forms()).toEqual([
      ['A ×2', 'residues'],
      ['Fc glycan', 'residues'],
      ['ATP', 'whole'],
    ]);
    const section = w.spec.root.value.children!.find((n) => n.id === 'views')!.children![0].children![0];
    expect(section.children!.filter((n) => n.kind === 'block').map((n) => [n.id, n.params])).toEqual([
      ['g', { form: 'residues' }],
    ]);
    // a unit the complex gains after the override: listed, derived
    const complex = w.spec.root.value.children!.find((n) => n.id === 'data')!.children!.find((n) => n.id === 'c')!;
    const a = complex.children!.find((n) => n.id === 'A')!;
    await w.spec.patch('Add a unit', {
      ops: [{ op: 'insert', parent: ['data', 'c'], place: { last: true }, node: { ...a, id: 'C' } }],
    });
    await w.settle();
    expect(forms()).toEqual([
      ['A ×3', 'residues'],
      ['Fc glycan', 'residues'],
      ['ATP', 'whole'],
    ]);
    // one undo step each: the unit, then the override
    const stack = (w.spec.history as UiHistory).stack;
    await stack.undo();
    await stack.undo();
    await w.settle();
    expect(forms()).toEqual([
      ['A ×2', 'residues'],
      ['Fc glycan', 'whole'],
      ['ATP', 'whole'],
    ]);
    await w.execute('block.form', { address: '/views/v/complex/l', form: 'atoms' });
    expect(w.view.resolved.value!.sections[0].blocks![2].axis.system).toBe('/data/c/l/(atoms)');
  });
});
