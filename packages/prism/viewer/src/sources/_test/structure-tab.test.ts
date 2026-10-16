/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { describe, expect, it } from 'vitest';
import { cbs, expectOk, headlessViewer, until } from '../../_test/harness.js';

// The Structure tab's Preset menu (state.md §1): the preset in effect is marked until the scene changes some other way.

describe('the presets model', () => {
  it('marks the preset applied, and not once Undo or a style has changed the scene', async () => {
    const v = await headlessViewer({ files: [cbs()] });
    const presets = v.views.presets;
    const idle = () => until(() => v.api.running.value.length === 0);
    expect(presets.current.value).toBeUndefined();

    presets.apply('illustrative');
    await until(() => presets.current.value === 'illustrative');
    await idle();
    // Nothing else changed: it stays.
    expect(presets.current.value).toBe('illustrative');

    expectOk(await v.api.history.undo());
    await until(() => presets.current.value === undefined);

    presets.apply('polymer-cartoon');
    await until(() => presets.current.value === 'polymer-cartoon');
    v.views.quickStyles.applyRepresentation(v.views.quickStyles.representations[0].id);
    await until(() => presets.current.value === undefined);
  });
});
