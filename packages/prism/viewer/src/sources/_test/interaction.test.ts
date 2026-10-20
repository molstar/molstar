/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { StructureElement } from '@molstar/model/model/structure';
import { describe, expect, it } from 'vitest';
import { cbs, crambin, expectOk, headlessViewer } from '../../_test/harness.js';

// What each structure has picked (Superpose reads it): follows the selection and the structures, and does not count a
// selection again for a change elsewhere in the tree.

describe('picked', () => {
  it('follows the selection; a change elsewhere keeps each row; a renamed structure gets a new row', async () => {
    const v = await headlessViewer({ files: [cbs()] });
    const picked = v.sources.interaction.picked;
    const seen: (typeof picked.value)[] = [];
    const sub = picked.subscribe((p) => seen.push(p));
    const data = () => v.plugin.managers.structure.hierarchy.current.structures[0].cell.obj!.data;
    v.plugin.managers.structure.selection.fromLoci('set', StructureElement.Loci.all(data()));
    const first = picked.value;
    expect(first).toHaveLength(1);
    expect(first[0].atoms).toBe(StructureElement.Loci.size(StructureElement.Loci.all(data())));

    expectOk(await v.api.styles.apply({ style: 'spacefill' }));
    expect(picked.value[0]).toBe(first[0]);

    // Another structure of the same name: the first one's label changes, so its row does.
    expectOk(await v.api.load.files({ files: [cbs()] }));
    expect(picked.value).toHaveLength(1);
    expect(picked.value[0].label).not.toBe(first[0].label);

    expectOk(await v.api.load.files({ files: [crambin()] }));
    v.plugin.managers.structure.selection.clear();
    expect(picked.value).toEqual([]);
    sub.unsubscribe();
  });
});
