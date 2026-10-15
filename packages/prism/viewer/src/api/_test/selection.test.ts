/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { StructureElement } from '@molstar/model/model/structure';
import { describe, expect, it } from 'vitest';
import { cbs, expectOk, headlessViewer } from '../../_test/harness.js';
import { createTargetResolver } from '../../sources/targets.js';
import type { PrismViewer } from '../../viewer.js';
import type { RunRecord } from '../run.js';

// `select`, `focus`, `highlight` and `camera` (operations.md §2, plan 2.3) on a headless viewer (no canvas) with 1CBS.

const selection = (v: PrismViewer) => v.plugin.managers.structure.selection;
const count = (v: PrismViewer) => selection(v).elementCount();
const steps = (v: PrismViewer) => v.api.history.steps.value.length;

async function setup() {
  return headlessViewer({ files: [cbs()] });
}

describe('select', () => {
  it('sets, adds, removes and intersects, with no history step', async () => {
    const v = await setup();
    const before = steps(v);
    const set = await v.api.select.set({ target: 'A:1-10' });
    expectOk(set);
    expect(set.did).toMatch(/^Selected 10 residues.* in chain A\.$/);
    expect(set.changed).toEqual({ selection: true });
    const ten = count(v);
    expect(ten).toBeGreaterThan(10);

    expectOk(await v.api.select.set({ target: 'A:5-20', mode: 'add' }));
    const twenty = count(v);
    expect(twenty).toBeGreaterThan(ten);

    expectOk(await v.api.select.set({ target: 'A:1-4', mode: 'remove' }));
    const removed = count(v);
    expect(removed).toBeLessThan(twenty);

    expectOk(await v.api.select.set({ target: 'A:10-12', mode: 'intersect' }));
    const three = count(v);
    expect(three).toBeLessThan(removed);
    expect(three).toBeGreaterThan(0);

    expectOk(await v.api.select.set({ target: 'ligand' }));
    expect(count(v)).not.toBe(three);
    expect(steps(v)).toBe(before);
  });

  it('says when it is already so, and when the target selects nothing', async () => {
    const v = await setup();
    expectOk(await v.api.select.set({ target: 'A:1-10' }));
    const again = await v.api.select.set({ target: 'A:1-10' });
    expectOk(again);
    expect(again.unchanged).toBe(true);
    const absent = await v.api.select.set({ target: 'A:50-60', mode: 'remove' });
    expectOk(absent);
    expect(absent.unchanged).toBe(true);
    const none = await v.api.select.set({ target: 'chain Z' });
    expect(none).toMatchObject({ ok: false, code: 'not-found' });
    expect(count(v)).toBeGreaterThan(0);
  });

  it('describes the selection as text and as a value', async () => {
    const v = await setup();
    const empty = await v.api.select.get();
    expectOk(empty);
    expect(empty.did).toBe('Nothing is selected.');
    expectOk(await v.api.select.set({ target: 'A:1-10' }));
    const got = await v.api.select.get();
    expectOk(got);
    expect(got.did).toMatch(/^Selected: 10 residues.* in chain A\.$/);
    expect(got.value).toHaveLength(1);
    expect(got.value?.[0].selected.residues).toBe(10);
  });

  it('selects all, clears, and inverts', async () => {
    const v = await setup();
    const all = await v.api.select.all();
    expectOk(all);
    const total = count(v);
    expect(total).toBeGreaterThan(1000);
    const again = await v.api.select.all();
    expectOk(again);
    expect(again.unchanged).toBe(true);

    const cleared = await v.api.select.clear();
    expectOk(cleared);
    expect(count(v)).toBe(0);
    const nothing = await v.api.select.clear();
    expectOk(nothing);
    expect(nothing.unchanged).toBe(true);

    expectOk(await v.api.select.set({ target: 'ligand' }));
    const ligand = count(v);
    expectOk(await v.api.select.invert());
    expect(count(v)).toBe(total - ligand);
    expectOk(await v.api.select.invert());
    expect(count(v)).toBe(ligand);
  });

  it('grows the selection around and expands it to whole residues and chains', async () => {
    const v = await setup();
    const none = await v.api.select.around({ radius: 5 });
    expect(none).toMatchObject({ ok: false, code: 'blocked' });

    expectOk(await v.api.select.set({ target: 'ligand' }));
    const ligand = count(v);
    const around = await v.api.select.around({ radius: 5 });
    expectOk(around);
    const grown = count(v);
    expect(grown).toBeGreaterThan(ligand);
    const atoms = await v.api.select.around({ radius: 8, asWholeResidues: false });
    expectOk(atoms);
    expect(count(v)).toBeGreaterThan(grown);

    expectOk(await v.api.select.set({ target: 'A:1@CA' }));
    expect(count(v)).toBe(1);
    expectOk(await v.api.select.expand({ to: 'residues' }));
    const residue = count(v);
    expect(residue).toBeGreaterThan(1);
    const same = await v.api.select.expand({ to: 'residues' });
    expectOk(same);
    expect(same.unchanged).toBe(true);
    expectOk(await v.api.select.expand({ to: 'chains' }));
    expect(count(v)).toBeGreaterThan(residue * 10);
    const got = await v.api.select.get();
    expectOk(got);
    expect(got.value?.[0].selected.chains).toEqual(['A']);
  });

  it('sets the selection mode and the picking level, and says when they are already so', async () => {
    const v = await setup();
    expectOk(await v.api.select.mode({ on: true }));
    expect(v.plugin.selectionMode).toBe(true);
    const mode = await v.api.select.mode({ on: true });
    expectOk(mode);
    expect(mode.unchanged).toBe(true);
    expectOk(await v.api.select.mode({ on: false }));
    expect(v.plugin.selectionMode).toBe(false);

    expectOk(await v.api.select.granularity({ granularity: 'chain' }));
    expect(v.plugin.managers.interactivity.props.granularity).toBe('chain');
    const level = await v.api.select.granularity({ granularity: 'chain' });
    expectOk(level);
    expect(level.unchanged).toBe(true);
    const wrong = await v.api.select.granularity({ granularity: 'nonsense' as never });
    expect(wrong).toMatchObject({ ok: false, code: 'invalid' });
    if (!wrong.ok) expect(wrong.hint).toMatch(/residue/);
  });

  it('takes loci from code in the page, directly', async () => {
    const v = await setup();
    const resolved = await createTargetResolver(v.plugin).resolve('A:20-24');
    if (!resolved.ok) throw new Error(resolved.message);
    const loci = resolved.parts[0].loci;
    const r = await v.api.select.set({ target: { loci } });
    expectOk(r);
    expect(r.did).toMatch(/^Selected 5 residues/);
    expect(count(v)).toBe(StructureElement.Loci.size(loci));
  });

  it('runs as the sequence link and is recorded with that origin', async () => {
    const v = await setup();
    const records: RunRecord[] = [];
    const sub = v.api.events.finished.subscribe((r) => records.push(r));
    const before = steps(v);
    const resolved = await createTargetResolver(v.plugin).resolve('A:20-24');
    if (!resolved.ok) throw new Error(resolved.message);
    const link = v.api.as({ kind: 'link' });
    expectOk(await link.select.set({ target: { loci: resolved.parts[0].loci } }));
    expectOk(await link.highlight.set({ target: { loci: resolved.parts[0].loci } }));
    expectOk(await link.highlight.clear());
    sub.unsubscribe();
    expect(records.map((r) => [r.operation, r.origin.kind])).toEqual([
      ['select.set', 'link'],
      ['highlight.set', 'link'],
      ['highlight.clear', 'link'],
    ]);
    expect(steps(v)).toBe(before);
  });
});

describe('focus', () => {
  it('sets and clears the focus without a canvas', async () => {
    const v = await setup();
    const focus = v.plugin.managers.structure.focus;
    expect(focus.current).toBeUndefined();
    const r = await v.api.focus.set({ target: 'ligand' });
    expectOk(r);
    expect(r.did).toMatch(/^Focused on /);
    expect(focus.current).toBeDefined();

    const got = await v.api.focus.get();
    expectOk(got);
    expect(got.value?.selected.residues).toBe(1);
    const again = await v.api.focus.set({ target: 'ligand' });
    expectOk(again);
    expect(again.unchanged).toBe(true);

    expectOk(await v.api.focus.clear());
    expect(focus.current).toBeUndefined();
    const nothing = await v.api.focus.clear();
    expectOk(nothing);
    expect(nothing.unchanged).toBe(true);
    const none = await v.api.focus.set({ target: 'chain Z' });
    expect(none).toMatchObject({ ok: false, code: 'not-found' });
  });
});

describe('highlight', () => {
  it('highlights and clears without a canvas', async () => {
    const v = await setup();
    const r = await v.api.highlight.set({ target: 'ligand' });
    expectOk(r);
    expect(r.changed).toEqual({ view: true });
    expectOk(await v.api.highlight.clear());
    expect(await v.api.highlight.set({ target: 'chain Z' })).toMatchObject({ ok: false, code: 'not-found' });
  });

  it("lights the sequence too when asked (the marks' hover), and clears it", async () => {
    const v = await setup();
    expectOk(await v.api.highlight.set({ target: 'ligand' }));
    expect(v.sources.marks.hover.value).toBeUndefined();
    expectOk(await v.api.highlight.set({ target: 'A:10', sequence: true }));
    expect(v.sources.marks.hover.value).toMatchObject({ owner: 'highlight' });
    expect(StructureElement.Loci.isEmpty(v.sources.marks.hover.value!.loci)).toBe(false);
    expectOk(await v.api.highlight.clear());
    expect(v.sources.marks.hover.value).toBeUndefined();
  });

  it("marks an owner's ranges in their tones on the representations, as view state, and takes them off", async () => {
    const v = await setup();
    const before = steps(v);
    const reprs = () =>
      v.plugin.managers.structure.hierarchy.current.structures[0].components.flatMap((c) =>
        c.representations.map((r) => r.cell.obj!.data.repr),
      );
    const layers = () => reprs().map((r) => r.state.overpaint.layers.length);
    const frame = () => new Promise<void>((r) => setTimeout(r, 20));
    expect(layers().every((n) => n === 0)).toBe(true);

    const r = await v.api.highlight.mark({
      owner: 'plot',
      marks: [
        { id: 'a', tone: 'highlight', label: 'Aligned A 1–10', target: 'A:1-10' },
        { id: 'b', tone: 'accent', target: 'A:50-60' },
      ],
    });
    expectOk(r);
    expect(r.did).toBe('Marked Aligned A 1–10; A:50-60.');
    expect(v.sources.marks.marks.value.map((m) => [m.owner, m.id, m.tone])).toEqual([
      ['plot', 'a', 'highlight'],
      ['plot', 'b', 'accent'],
    ]);
    await frame();
    // The polymer's representation has both; the ligand's neither.
    expect(Math.max(...layers())).toBe(2);
    expect(steps(v)).toBe(before);

    // Another owner's marks stay when one clears its own.
    expectOk(await v.api.highlight.mark({ owner: 'other', marks: [{ id: 'x', tone: 'accent', target: 'ligand' }] }));
    expectOk(await v.api.highlight.unmark({ owner: 'plot' }));
    expect(v.sources.marks.marks.value.map((m) => m.owner)).toEqual(['other']);
    expectOk(await v.api.highlight.unmark({ owner: 'other' }));
    await frame();
    expect(layers().every((n) => n === 0)).toBe(true);

    // The rest greyed: a layer under the marks on every representation, the ligand's (none of the marks) too.
    expectOk(
      await v.api.highlight.mark({
        owner: 'plot',
        dimRest: true,
        marks: [{ id: 'a', tone: 'accent', target: 'A:1-10' }],
      }),
    );
    await frame();
    expect(Math.max(...layers())).toBe(2);
    expect(Math.min(...layers())).toBe(1);
    expectOk(await v.api.highlight.unmark({ owner: 'plot' }));
    await frame();
    expect(layers().every((n) => n === 0)).toBe(true);

    expect(await v.api.highlight.mark({ owner: '', marks: [] })).toMatchObject({ ok: false, code: 'invalid' });
    expect(
      await v.api.highlight.mark({ owner: 'plot', marks: [{ id: 'z', tone: 'accent', target: 'chain Z' }] }),
    ).toMatchObject({ ok: false, code: 'not-found' });
  });
});

describe('camera', () => {
  it('needs a canvas to frame, zoom and orient, and describes itself without one', async () => {
    const v = await setup();
    expect(await v.api.camera.frame({ target: 'ligand' })).toMatchObject({ ok: false, code: 'unavailable' });
    expect(await v.api.camera.zoom({ factor: 2 })).toMatchObject({ ok: false, code: 'unavailable' });
    expect(await v.api.camera.orient()).toMatchObject({ ok: false, code: 'unavailable' });
    const got = await v.api.camera.get();
    expectOk(got);
    expect(got.did).toMatch(/No canvas/);
  });
});

describe('the selection history (gate G5)', () => {
  it('gets what was added, not the whole selection', async () => {
    const v = await headlessViewer({ files: [cbs()] });
    expectOk(await v.api.select.set({ target: 'ligand' }));
    expectOk(await v.api.select.set({ target: 'A:58', mode: 'add' }));
    const history = v.plugin.managers.structure.selection.additionsHistory;
    const latest = history[0].loci;
    expect(StructureElement.Loci.size(latest)).toBeLessThan(v.plugin.managers.structure.selection.elementCount());
    expect(history).toHaveLength(2);
  });
});
