/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { StructureElement } from '@molstar/model/model/structure';
import { describe, expect, it } from 'vitest';
import { cbs, crambin, expectOk, headlessViewer, until } from '../../_test/harness.js';
import type { RunRecord } from '../../api/run.js';

// The sequence on the API (sequence.md §2–§7, plan 2.4), headless.

describe('the sequence', () => {
  it('follows the loaded structure, and shows its document', async () => {
    const v = await headlessViewer({ files: [cbs()] });
    const seq = v.views.sequence!;
    expect(seq.source.structure.value?.title).toBe('1CBS');
    const viewer = seq.viewers.get();
    await until(() => !!viewer.structure.value);
    expect(viewer.structure.value?.title).toBe('1CBS');
  });

  it('follows another structure when asked, and the rule again', async () => {
    const v = await headlessViewer({ files: [cbs(), crambin()] });
    const seq = v.views.sequence!;
    expect(seq.source.structure.value?.title).toBe('1CBS');
    expectOk(await v.api.sequence!.follow({ structure: '1crn' }));
    expect(seq.source.structure.value?.title).toBe('1CRN');
    expectOk(await v.api.sequence!.follow({ structure: '' }));
    expect(seq.source.structure.value?.title).toBe('1CBS');
  });

  it('changes its style as a step, and Undo brings it back (its state is in the checkpoint)', async () => {
    const v = await headlessViewer({ files: [cbs()] });
    const viewer = v.views.sequence!.viewers.get();
    expect(viewer.style.value).toBe('tracks');
    expectOk(await v.api.sequence!.style({ style: 'text' }));
    expect(viewer.style.value).toBe('compact');
    const steps = v.managers.history.steps.value.length;
    expect(v.managers.history.steps.value[0].label).toMatch(/text/);
    // Already so: nothing changes, no step.
    expect(await v.api.sequence!.style({ style: 'text' })).toMatchObject({ ok: true, unchanged: true });
    expect(v.managers.history.steps.value).toHaveLength(steps);
    await v.api.history.undo();
    expect(viewer.style.value).toBe('tracks');
  });

  it("writes into Mol* only through operations, with origin 'link'", async () => {
    const v = await headlessViewer({ files: [cbs()] });
    const seq = v.views.sequence!;
    const records: RunRecord[] = [];
    v.api.events.finished.subscribe((r) => records.push(r));
    const structure = seq.source.structure.value!.structure;
    const view = seq.source.view3d(structure);
    const loci = StructureElement.Loci.firstResidue(StructureElement.Loci.all(structure));
    view.select(loci);
    await until(() => records.some((r) => r.operation === 'select.set'));
    view.highlight(loci);
    // Sent back to back: they apply in order (a selection then its clearing never ends selected).
    view.select(loci);
    view.select(null);
    await until(() => records.at(-1)?.operation === 'select.clear');
    const ops = records.map((r) => `${r.origin.kind}:${r.operation}`);
    expect(ops[0]).toBe('link:select.set');
    expect(ops).toContain('link:highlight.set');
    expect(ops.at(-1)).toBe('link:select.clear');
    expect(v.plugin.managers.structure.selection.elementCount()).toBe(0);
    // View changes: no history steps.
    expect(v.managers.history.steps.value).toHaveLength(1);
  });

  it("colors the 3D view from a track's menu as a step", async () => {
    const v = await headlessViewer({ files: [cbs()] });
    const seq = v.views.sequence!;
    await seq.source.actions.colorStructure('secondary-structure');
    expect(seq.source.actions.coloring?.()).toBe('secondary-structure');
    expect(v.managers.history.steps.value[0].origin.kind).toBe('link');
  });

  it('is absent when the spec has none', async () => {
    const v = await headlessViewer({ spec: { sequence: false } });
    expect(v.views.sequence).toBeUndefined();
    expect(v.api.sequence).toBeUndefined();
  });
});
