/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { describe, expect, it } from 'vitest';
import { cbs, expectOk, headlessViewer, representationTypes } from '../../_test/harness.js';

// `styles` (api.md §2): quick styles.

describe('styles', () => {
  it('applies representation styles as one step each, and looks', async () => {
    const v = await headlessViewer({ files: [cbs()] });
    const before = representationTypes(v);
    expectOk(await v.api.styles.apply({ style: 'surface' }));
    expect(representationTypes(v)).toContain('molecular-surface');
    expect(v.managers.history.steps.value).toHaveLength(2);

    expectOk(await v.api.styles.apply({ style: 'illustrative' }));
    expect(v.plugin.managers.structure.component.state.options.ignoreLight).toBe(true);
    const again = await v.api.styles.apply({ style: 'illustrative' });
    expect(again).toMatchObject({ ok: true, unchanged: true });

    // A representation style keeps the look in effect.
    expectOk(await v.api.styles.apply({ style: 'default' }));
    expect(v.plugin.managers.structure.component.state.options.ignoreLight).toBe(true);
    expect(representationTypes(v)).toEqual(before);

    expectOk(await v.api.styles.apply({ style: 'plain' }));
    expect(v.plugin.managers.structure.component.state.options.ignoreLight).toBe(false);
  });

  it('keeps a plain look on the representations of a preset that brings its own lighting', async () => {
    const v = await headlessViewer({ files: [cbs()] });
    expectOk(await v.api.styles.apply({ style: 'spacefill' }));
    // the spacefill preset ignores light; the look in effect (plain) is on every representation, not an override
    const reprs = v.plugin.managers.structure.hierarchy.current.structures[0].components.flatMap(
      (c) => c.representations,
    );
    expect(reprs.length).toBeGreaterThan(0);
    for (const r of reprs) expect((r.cell.transform.params as any).type.params.ignoreLight).toBe(false);
    expect(await v.api.look.get({})).toMatchObject({ value: { overrides: [] } });
  });

  it('refuses an unknown style with the list', async () => {
    const v = await headlessViewer({ files: [cbs()] });
    const r = await v.api.run('styles.apply', { style: 'neon' });
    expect(r).toMatchObject({ ok: false, code: 'invalid' });
    if (!r.ok) expect(r.hint).toMatch(/cartoon/);
  });
});
