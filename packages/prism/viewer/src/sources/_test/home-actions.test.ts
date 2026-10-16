/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { StateAction } from '@molstar/core/state';
import { PluginStateObject } from '@molstar/plugin/state/objects';
import { describe, expect, it } from 'vitest';
import { HomeFeature } from '../../features/home/index.js';
import { headlessViewer } from '../../_test/harness.js';

// Home's actions follow the registry, and an action registered later leaves the forms of the others as they are.

describe('Home actions', () => {
  it('keep their forms when another action is registered or removed', async () => {
    const v = await headlessViewer({ extensions: [HomeFeature] });
    const before = v.views.homeActions.value;
    expect(before.length).toBeGreaterThan(0);

    const Late = StateAction.build({
      from: PluginStateObject.Root,
      display: { name: 'Late action' },
    })(async () => {});
    v.plugin.state.data.actions.add(Late);
    const during = v.views.homeActions.value;
    expect(during.map((a) => a.title)).toContain('Late action');
    for (const a of before) expect(during.find((d) => d.id === a.id)).toBe(a);

    v.plugin.state.data.actions.remove(Late);
    const after = v.views.homeActions.value;
    expect(after.map((a) => a.title)).not.toContain('Late action');
    expect(after).toHaveLength(before.length);
    for (const a of before) expect(after.find((d) => d.id === a.id)).toBe(a);
  });
});
