/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { EmptyLoci } from '@molstar/model/model/loci';
import { StructureElement } from '@molstar/model/model/structure';
import { describe, expect, it } from 'vitest';
import { cbs, crambin, expectOk, headlessViewer, until } from '../../_test/harness.js';
import type { RunRecord } from '../../api/run.js';

// The sequence on the API (sequence.md §2–§7), headless.

describe('the sequence', () => {
  it('follows the loaded structure, and shows its document', async () => {
    const v = await headlessViewer({ files: [cbs()] });
    const seq = v.views.sequence!;
    expect(seq.source.structure.value?.title).toBe('1CBS');
    const viewer = (await seq.load()).viewers.get();
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
    const viewer = (await v.views.sequence!.load()).viewers.get();
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

  it('reads the coloring of the structure it shows, and colors that one only', async () => {
    const v = await headlessViewer({ files: [cbs(), crambin()] });
    const seq = v.views.sequence!;
    expect(seq.source.structure.value?.title).toBe('1CBS');
    await seq.source.actions.colorStructure('secondary-structure');
    // The other structure is colored the way it was: the menu still reads the shown one's.
    expect(seq.source.actions.coloring?.()).toBe('secondary-structure');
    expectOk(await v.api.sequence!.follow({ structure: '1crn' }));
    expect(seq.source.actions.coloring?.()).not.toBe('secondary-structure');
  });

  it('follows the structure whose selection grew, not the last one with a selection', async () => {
    const v = await headlessViewer({ files: [cbs(), crambin()] });
    const seq = v.views.sequence!;
    expectOk(await v.api.select.set({ target: 'A:1-5', structure: '1crn' }));
    expect(seq.source.structure.value?.title).toBe('1CRN');
    expectOk(await v.api.select.set({ target: 'A:1-5', structure: '1cbs' }));
    expect(seq.source.structure.value?.title).toBe('1CBS');
    // Both selected: the one that grows is shown, whatever its place in the list.
    expectOk(await v.api.select.set({ target: 'A:10-12', structure: '1crn', mode: 'add' }));
    expect(seq.source.structure.value?.title).toBe('1CRN');
    expectOk(await v.api.select.set({ target: 'A:10-12', structure: '1cbs', mode: 'add' }));
    expect(seq.source.structure.value?.title).toBe('1CBS');
  });

  it('passes the 3D hover on once per element', async () => {
    const v = await headlessViewer({ files: [cbs()] });
    const seq = v.views.sequence!;
    const structure = seq.source.structure.value!.structure;
    const view = seq.source.view3d(structure);
    const first = StructureElement.Loci.firstResidue(StructureElement.Loci.all(structure));
    const seen: unknown[] = [];
    const stop = view.onHover((l) => seen.push(l));
    // The behavior replays its last value on subscribing.
    seen.length = 0;
    const hover = (loci: StructureElement.Loci | EmptyLoci) =>
      v.plugin.behaviors.interaction.hover.next({
        current: { loci },
        buttons: 0,
        button: 0,
        modifiers: { shift: false, control: false, alt: false, meta: false },
      } as never);
    hover(first);
    hover(StructureElement.Loci.remap(first, structure));
    hover(first);
    expect(seen).toHaveLength(1);
    hover(EmptyLoci);
    hover(first);
    expect(seen).toHaveLength(3);
    stop();
  });

  it('is absent when the spec has none', async () => {
    const v = await headlessViewer({ spec: { sequence: false } });
    expect(v.views.sequence).toBeUndefined();
    expect(v.api.sequence).toBeUndefined();
  });

  it('loads its runtime only when a structure is shown, and applies a saved state given before', async () => {
    // A state saved with the text style, a structure followed.
    const a = await headlessViewer({ files: [cbs()] });
    expectOk(await a.api.sequence!.style({ style: 'text' }));
    const snapshot = a.plugin.state.getSnapshot({ camera: false, canvas3d: false });

    const v = await headlessViewer();
    const seq = v.views.sequence!;
    expect(seq.viewers.value).toBeUndefined();
    // The state arrives before the runtime: kept, then applied.
    await v.plugin.state.setSnapshot(snapshot);
    await until(() => !!seq.viewers.value);
    const viewer = seq.viewers.value!.get();
    await until(() => viewer.style.value === 'compact');
    expect(
      v.plugin.state.getSnapshot({ camera: false, canvas3d: false }).custom?.['prism-viewer.sequence'],
    ).toMatchObject({ style: 'compact' });
  });
});
