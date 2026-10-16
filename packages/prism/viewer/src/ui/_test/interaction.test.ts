/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { describe, expect, it } from 'vitest';
import { headlessViewer } from '../../_test/harness.js';
import { selectionModel } from '../interaction.js';

// The selection bar's actions follow the viewer's features.

describe('the selection bar', () => {
  it('offers Component only where the viewer has the component area', async () => {
    const v = await headlessViewer();
    const person = v.api.as({ kind: 'person' });
    const ids = (api: typeof person) => selectionModel(v.sources.interaction, api).actions.map((a) => a.id);
    expect(ids(person)).toContain('component');
    const without = new Proxy(person, {
      get: (target, key) =>
        key === 'area' ? (id: string) => (id === 'component' ? undefined : target.area(id)) : (target as never)[key],
    });
    expect(ids(without)).not.toContain('component');
    expect(ids(without)).toContain('focus');
  });
});
