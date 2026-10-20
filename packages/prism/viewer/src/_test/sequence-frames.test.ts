/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { type Structure, StructureElement } from '@molstar/model/model/structure';
import { lociSize } from '@molstar/seq-viewer';
import { describe, expect, it, vi } from 'vitest';
import { ensemble } from '../api/_test/trajectory-fixtures.js';
import { linkSequence } from '../sequence/link.js';
import { cbs, expectOk, headlessViewer, until } from './harness.js';

// A trajectory's frames in the sequence view (sequence.md §2): a frame of the structure shown keeps the document and
// the link's systems; the link maps onto the frame drawn.

const drawn = (v: Awaited<ReturnType<typeof headlessViewer>>) =>
  v.plugin.managers.structure.hierarchy.current.structures[0].cell.obj!.data;

describe('the sequence over frames', () => {
  it('keeps the document and the link over frame steps; maps the selection onto the frame drawn', async () => {
    const v = await headlessViewer({ files: [ensemble()] });
    const seqView = v.views.sequence!;
    const viewer = (await seqView.load()).viewers.get();
    await until(() => !!viewer.structure.value && !!viewer.seq.graph.assembled.value);
    const shown = viewer.structure.value!;
    const set = viewer.seq.graph.assembled.value!.set;
    const document = viewer.seq.document.shown.value;
    const register = vi.spyOn(viewer.seq.systems, 'register');
    const views: Structure[] = [];
    const stop = linkSequence(viewer, (s) => {
      views.push(s);
      return seqView.source.view3d(s);
    });
    expect(views).toEqual([drawn(v)]);
    const registered = register.mock.calls.length;
    expect(registered).toBeGreaterThan(0);

    for (const by of [1, 1, -1]) {
      const before = drawn(v);
      expectOk(await v.api.trajectory.frame({ by }));
      expect(drawn(v)).not.toBe(before);
      await until(() => viewer.frame.value === drawn(v));
    }
    // No new document, no band emptied and refilled, no systems registered again.
    expect(viewer.structure.value).toBe(shown);
    expect(viewer.seq.graph.assembled.value!.set).toBe(set);
    expect(viewer.seq.document.shown.value).toBe(document);
    expect(register.mock.calls.length).toBe(registered);
    // The views follow the frame drawn: a 3D selection on it maps into the sequence.
    expect(views.at(-1)).toBe(drawn(v));
    v.plugin.managers.structure.selection.fromLoci('set', StructureElement.Loci.all(drawn(v)));
    await until(() => (lociSize(viewer.seq.selection.current.value.seq)?.length ?? 0) > 0);
    stop();
  });

  it('builds a new document for another structure', async () => {
    const v = await headlessViewer({ files: [ensemble()] });
    const seqView = v.views.sequence!;
    const viewer = (await seqView.load()).viewers.get();
    await until(() => !!viewer.structure.value);
    const first = viewer.structure.value;
    expectOk(await v.api.load.files({ files: [cbs()] }));
    const other = v.plugin.managers.structure.hierarchy.current.structures[1].cell.transform.ref;
    seqView.source.setFollowed(other);
    await until(() => !!viewer.structure.value && viewer.structure.value !== first);
    expect(viewer.frame.value).toBe(viewer.structure.value!.structure);
    expect(StructureElement.Loci.is(StructureElement.Loci.all(viewer.frame.value!))).toBe(true);
  });
});
