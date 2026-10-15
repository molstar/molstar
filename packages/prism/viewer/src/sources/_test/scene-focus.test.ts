/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { describe, expect, it } from 'vitest';
import { cbs, headlessViewer, until } from '../../_test/harness.js';

// The focus group: always under a structure (the tree does not shift as the focus comes and goes), holding the focus
// behavior's components while something is focused.

describe('the focus group', () => {
  it('is there before, during and after a focus, with the target and the surroundings in it', async () => {
    const v = await headlessViewer({ files: [cbs()] });
    const focus = () => v.views.scene.roots.value[0].children.value.find((c) => c.kind === 'focus')!;
    expect(focus().description.value).toBe('none');
    expect(focus().children.value).toEqual([]);
    const components = () => v.views.scene.roots.value[0].children.value.filter((c) => c.kind === 'component');
    const before = components().length;

    const r = await v.api.focus.set({ target: 'ligand' });
    expect(r.ok).toBe(true);
    await until(() => focus().children.value.length === 2);
    expect(focus().children.value.map((c) => c.label.value)).toEqual(['Target', 'Surroundings (5 Å)']);
    expect(focus().description.value).toMatch(/\d+ atoms/);
    // The focus parts are not among the structure's own components.
    expect(components().length).toBe(before);

    await v.api.focus.clear({});
    await until(() => focus().children.value.length === 0);
    expect(focus().description.value).toBe('none');
  });
});
