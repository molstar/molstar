/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { describe, expect, it } from 'vitest';
import { cbs, crambin, expectOk, headlessViewer, representationTypes } from '../../_test/harness.js';
import type { PrismViewer } from '../../viewer.js';
import type { StructureComponents } from '../component.js';

// `component` (api.md §2) on a headless viewer with 1CBS.

const steps = (v: PrismViewer) => v.managers.history.steps.value.length;
const components = (v: PrismViewer) =>
  v.plugin.managers.structure.hierarchy.current.structures.flatMap((s) => s.components);
const named = (v: PrismViewer, label: string) =>
  components(v).filter((c) => (c.cell.params?.values as any).label === label);
const atoms = (v: PrismViewer, label: string) => named(v, label)[0]?.cell.obj?.data.elementCount;

async function withLigand(files = [cbs()]) {
  const v = await headlessViewer({ files });
  const r = await v.api.component.add({ target: 'ligand', label: 'Retinoic acid', representation: 'ball-and-stick' });
  expectOk(r);
  return { v, r };
}

describe('component.add', () => {
  it('makes a component of a target with a representation, one undoable step', async () => {
    const v = await headlessViewer({ files: [cbs()] });
    const before = representationTypes(v);
    const stepsBefore = steps(v);
    const r = await v.api.component.add({ target: 'ligand', label: 'Retinoic acid', representation: 'ball-and-stick' });
    expectOk(r);
    expect(r.changed.nodes).toHaveLength(1);
    expect(named(v, 'Retinoic acid')).toHaveLength(1);
    const c = named(v, 'Retinoic acid')[0];
    expect(c.cell.transform.ref).toBe(r.changed.nodes![0]);
    expect(c.representations.map((x) => (x.cell.transform.params as any).type.name)).toEqual(['ball-and-stick']);
    expect(representationTypes(v)).toHaveLength(before.length + 1);
    expect(steps(v)).toBe(stepsBefore + 1);

    await v.api.history.undo();
    expect(named(v, 'Retinoic acid')).toHaveLength(0);
    expect(representationTypes(v)).toEqual(before);
  });

  it('names the component after the target when there is no label, and uses ball-and-stick by default', async () => {
    const v = await headlessViewer({ files: [cbs()] });
    const r = await v.api.component.add({ target: 'resn REA' });
    expectOk(r);
    const c = named(v, 'resn REA')[0];
    expect(c).toBeDefined();
    expect(c.representations.map((x) => (x.cell.transform.params as any).type.name)).toEqual(['ball-and-stick']);
  });

  it('says it already exists, without a step, and makes another on copy', async () => {
    const { v } = await withLigand();
    const n = components(v).length;
    const stepsBefore = steps(v);
    // the same target, whatever its name: it is the same component
    const again = await v.api.component.add({ target: 'ligand', label: 'Something else' });
    expect(again).toMatchObject({ ok: true, unchanged: true });
    expect(components(v)).toHaveLength(n);
    expect(steps(v)).toBe(stepsBefore);
    // the same name, another target: another component (decision 30: identity is the structure and the target)
    expectOk(await v.api.component.add({ target: 'resi 10-20', label: 'Retinoic acid' }));
    expect(components(v)).toHaveLength(n + 1);

    const copy = await v.api.component.add({ target: 'ligand', label: 'Another' }, { copy: true });
    expectOk(copy);
    expect(components(v)).toHaveLength(n + 2);
    expect(steps(v)).toBe(stepsBefore + 2);
  });

  it('adds a component with no representation on none', async () => {
    const v = await headlessViewer({ files: [cbs()] });
    const before = representationTypes(v);
    expectOk(await v.api.component.add({ target: 'ligand', label: 'Bare', representation: 'none' }));
    expect(named(v, 'Bare')[0].representations).toHaveLength(0);
    expect(representationTypes(v)).toEqual(before);
  });

  it('fails when the target selects nothing, with nothing changed', async () => {
    const v = await headlessViewer({ files: [cbs()] });
    const n = components(v).length;
    const stepsBefore = steps(v);
    const r = await v.api.component.add({ target: 'resn ZZZ' });
    expect(r).toMatchObject({ ok: false, code: 'not-found' });
    if (!r.ok) expect(r.message).toMatch(/selects nothing/);
    expect(components(v)).toHaveLength(n);
    expect(steps(v)).toBe(stepsBefore);
  });

  it('fails on a representation that does not exist, naming the valid ones', async () => {
    const v = await headlessViewer({ files: [cbs()] });
    const r = await v.api.component.add({ target: 'ligand', representation: 'blobs' });
    expect(r).toMatchObject({ ok: false, code: 'invalid' });
    if (!r.ok) expect(r.hint).toMatch(/ball-and-stick/);
  });

  it('fails on a representation that does not apply to what is selected', async () => {
    const v = await headlessViewer({ files: [cbs()] });
    const r = await v.api.component.add({ target: 'ligand', representation: 'cartoon' });
    expect(r).toMatchObject({ ok: false, code: 'invalid' });
    if (!r.ok) expect(r.hint).toMatch(/ball-and-stick/);
  });

  it('asks which structure for a component name that is in both, and names the structure otherwise', async () => {
    const v = await headlessViewer({ files: [cbs(), crambin()] });
    // adds to every structure the target selects something in; a component name found in both is ambiguous
    const r = await v.api.component.rename({ component: 'Polymer', label: 'Chains' });
    expect(r).toMatchObject({ ok: false, code: 'ambiguous' });
    if (!r.ok) expect(r.hint).toMatch(/1CBS/i);
    expectOk(await v.api.component.rename({ component: 'Polymer', label: 'Chains', structure: '1cbs' }));
  });
});

describe('component.rename, visibility and remove', () => {
  it('renames a component, and undoes; the same name is already so', async () => {
    const { v } = await withLigand();
    const stepsBefore = steps(v);
    expectOk(await v.api.component.rename({ component: 'Retinoic acid', label: 'RA' }));
    expect(named(v, 'RA')).toHaveLength(1);
    expect(named(v, 'Retinoic acid')).toHaveLength(0);
    expect(steps(v)).toBe(stepsBefore + 1);

    expect(await v.api.component.rename({ component: 'RA', label: 'RA' })).toMatchObject({ ok: true, unchanged: true });
    expect(steps(v)).toBe(stepsBefore + 1);

    await v.api.history.undo();
    expect(named(v, 'Retinoic acid')).toHaveLength(1);
    expect(named(v, 'RA')).toHaveLength(0);
  });

  it('asks which component when none is named and the structure has several', async () => {
    const { v } = await withLigand();
    expect(components(v).length).toBeGreaterThan(1);
    const r = await v.api.component.rename({ component: '', label: 'X' });
    expect(r).toMatchObject({ ok: false, code: 'ambiguous' });
    if (!r.ok) expect(r.hint).toMatch(/Retinoic acid \(/);
    // The same rule in the representation area, where an empty component used to be `invalid`.
    expect(await v.api.representation.remove({ component: '' })).toMatchObject({ ok: false, code: 'ambiguous' });
  });

  it('finds a component by state ref or case-insensitive name, and lists the components when it does not', async () => {
    const { v, r } = await withLigand();
    expectOk(await v.api.component.rename({ component: r.changed.nodes![0], label: 'RA' }));
    expectOk(await v.api.component.rename({ component: 'ra', label: 'Ret' }));
    const missing = await v.api.component.rename({ component: 'nope', label: 'x' });
    expect(missing).toMatchObject({ ok: false, code: 'not-found' });
    if (!missing.ok) expect(missing.hint).toMatch(/Ret/);
    expect(await v.api.component.rename({ component: 'Ret', label: '  ' })).toMatchObject({
      ok: false,
      code: 'invalid',
    });
  });

  it('hides and shows a component, only when it differs, and undoes', async () => {
    const { v } = await withLigand();
    const hidden = () => named(v, 'Retinoic acid')[0].cell.state.isHidden;
    const stepsBefore = steps(v);
    expect(await v.api.component.visibility({ component: 'Retinoic acid', visible: true })).toMatchObject({
      ok: true,
      unchanged: true,
    });
    expect(steps(v)).toBe(stepsBefore);

    expectOk(await v.api.component.visibility({ component: 'Retinoic acid', visible: false }));
    expect(hidden()).toBe(true);
    expect(steps(v)).toBe(stepsBefore + 1);
    expect(await v.api.component.visibility({ component: 'Retinoic acid', visible: false })).toMatchObject({
      ok: true,
      unchanged: true,
    });

    await v.api.history.undo();
    expect(hidden()).toBeFalsy();
  });

  it('removes a component with its representations, and undoes', async () => {
    const { v } = await withLigand();
    const withRa = representationTypes(v);
    const stepsBefore = steps(v);
    const r = await v.api.component.remove({ component: 'Retinoic acid' });
    expectOk(r);
    expect(named(v, 'Retinoic acid')).toHaveLength(0);
    expect(representationTypes(v)).toHaveLength(withRa.length - 1);
    expect(steps(v)).toBe(stepsBefore + 1);

    const gone = await v.api.component.remove({ component: 'Retinoic acid' });
    expect(gone).toMatchObject({ ok: false, code: 'not-found' });
    if (!gone.ok) expect(gone.hint).toMatch(/Polymer/);

    await v.api.history.undo();
    expect(named(v, 'Retinoic acid')).toHaveLength(1);
    expect(representationTypes(v)).toEqual(withRa);
  });
});

describe('component.modify', () => {
  it('takes atoms out of a component, adds them back, and undoes', async () => {
    const v = await headlessViewer({ files: [cbs()] });
    expectOk(await v.api.component.add({ target: 'protein', label: 'Chain', representation: 'ball-and-stick' }));
    const all = atoms(v, 'Chain')!;
    const stepsBefore = steps(v);

    const sub = await v.api.component.modify({ component: 'Chain', target: 'resi 1-10', mode: 'remove' });
    expectOk(sub);
    const less = atoms(v, 'Chain')!;
    expect(less).toBeLessThan(all);
    expect(steps(v)).toBe(stepsBefore + 1);

    // taking out what is not there changes nothing
    expect(await v.api.component.modify({ component: 'Chain', target: 'resi 1-10', mode: 'remove' })).toMatchObject({
      ok: true,
      unchanged: true,
    });
    expect(steps(v)).toBe(stepsBefore + 1);

    expectOk(await v.api.component.modify({ component: 'Chain', target: 'resi 1-10', mode: 'add' }));
    expect(atoms(v, 'Chain')).toBe(all);
    expect(await v.api.component.modify({ component: 'Chain', target: 'resi 1-10', mode: 'add' })).toMatchObject({
      ok: true,
      unchanged: true,
    });

    await v.api.history.undo();
    expect(atoms(v, 'Chain')).toBe(less);
  });

  it('keeps only the overlap on intersect, and refuses an overlap of nothing', async () => {
    const v = await headlessViewer({ files: [cbs()] });
    expectOk(await v.api.component.add({ target: 'protein', label: 'Chain', representation: 'ball-and-stick' }));
    expectOk(await v.api.component.modify({ component: 'Chain', target: 'resi 1-10', mode: 'intersect' }));
    const some = atoms(v, 'Chain')!;
    expect(some).toBeGreaterThan(0);
    expect(some).toBeLessThan(200);
    const none = await v.api.component.modify({ component: 'Chain', target: 'resi 100-110', mode: 'intersect' });
    expect(none).toMatchObject({ ok: false, code: 'invalid' });
    expect(atoms(v, 'Chain')).toBe(some);
  });

  it('fails when the target selects nothing', async () => {
    const { v } = await withLigand();
    const r = await v.api.component.modify({ component: 'Retinoic acid', target: 'resn ZZZ', mode: 'add' });
    expect(r).toMatchObject({ ok: false, code: 'not-found' });
  });
});

describe('component.list', () => {
  it('reads the components of each structure with their representations', async () => {
    const { v, r } = await withLigand();
    const stepsBefore = steps(v);
    const list = await v.api.component.list({});
    expectOk(list);
    const value = list.value as StructureComponents[];
    expect(value).toHaveLength(1);
    expect(value[0].structure.label).toBe('1CBS');
    const ra = value[0].components.find((c) => c.label === 'Retinoic acid');
    expect(ra).toMatchObject({ ref: r.changed.nodes![0], visible: true });
    expect(ra?.representations.map((x) => x.type)).toEqual(['ball-and-stick']);
    expect(ra?.atoms).toBeGreaterThan(0);
    expect(list.did).toMatch(/Retinoic acid.*ball-and-stick/);
    expect(steps(v)).toBe(stepsBefore);
  });

  it('names the structure it was asked about, and fails on one that is not loaded', async () => {
    const v = await headlessViewer({ files: [cbs(), crambin()] });
    const one = await v.api.component.list({ structure: '1crn' });
    expectOk(one);
    expect((one.value as StructureComponents[]).map((s) => s.structure.label)).toEqual(['1CRN']);
    expect(await v.api.component.list({ structure: 'zzz' })).toMatchObject({ ok: false, code: 'not-found' });
  });
});
