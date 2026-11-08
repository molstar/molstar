/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Vec3 } from '@molstar/core/math/linear-algebra';
import { describe, expect, it } from 'vitest';
import { cbs, expectOk, headlessViewer } from '../../../_test/harness.js';
import { findCuts } from '../../../sources/clipping/cuts.js';
import type { PrismViewer } from '../../../viewer.js';
import { ExportSessionFeature } from '../../export-session/index.js';

// A session keeps the cuts: they are the state tree (clipping.md §1), and their drawings cells in it (§6).

/** The drawings: each cell tagged by the cut it draws (its representation below it has the shared tag only) */
const drawings = (v: PrismViewer) =>
  [...v.plugin.state.data.cells.values()].filter((c) => c.transform.tags?.some((t) => t.startsWith('clip-volume-')));

describe('a session with cuts', () => {
  it('opens with its cuts and their drawings, which stay editable', async () => {
    const v = await headlessViewer({ extensions: [ExportSessionFeature], files: [cbs()] });
    const added = await v.api.clipping.add({ shape: 'cube', target: 'polymer', keepWhole: true });
    expectOk(added);
    expectOk(await v.api.clipping.show({ cut: 1, on: true, options: { wallOpacity: 0.4 } }));
    expectOk(await v.api.clipping.update({ cut: 1, position: [5, 6, 7] }));
    const [saved] = findCuts(v.plugin);
    const session = await v.api.session.save!({ kind: 'session' });
    expectOk(session);

    const opened = await headlessViewer({ extensions: [ExportSessionFeature] });
    expectOk(await opened.api.session.open({ file: new File([session.value!.blob], 'cuts.molx') }));
    const cuts = findCuts(opened.plugin);
    expect(cuts).toHaveLength(1);
    expect(cuts[0].object.tag).toBe(saved.object.tag);
    expect(cuts[0].targets.map((t) => t.ref)).toEqual(saved.targets.map((t) => t.ref));
    expect(Vec3.exactEquals(cuts[0].object.position, Vec3.create(5, 6, 7))).toBe(true);
    expect(cuts[0].kept).toEqual(saved.kept);
    // one drawing, under its cut's representation, with its options, built
    expect(drawings(opened)).toHaveLength(1);
    const drawing = opened.plugin.state.data.cells.get(cuts[0].drawing!)!;
    expect(drawing.transform.parent).toBe(cuts[0].targets[0].ref);
    expect(drawing.transform.params).toMatchObject({ wallOpacity: 0.4 });
    expect(drawing.obj).toBeDefined();

    // the same cut as before: drawn once, moved with the drawing, removed with it
    expect(await opened.api.clipping.show({ cut: 1, on: true })).toMatchObject({ ok: true, unchanged: true });
    const before = drawing.obj;
    expectOk(await opened.api.clipping.update({ cut: 1, position: [1, 2, 3] }));
    expect(opened.plugin.state.data.cells.get(findCuts(opened.plugin)[0].drawing!)!.obj).not.toBe(before);
    expect(drawings(opened)).toHaveLength(1);
    expectOk(await opened.api.clipping.remove({ cut: 1 }));
    expect(drawings(opened)).toHaveLength(0);
  });
});
