/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { describe, expect, it } from 'vitest';
import { observableState } from '@molstar/prism-ui/core';
import type { SelectionAction } from '../../model/selection.js';
import { headlessViewer } from '../../_test/harness.js';
import { selectionModel } from '../interaction.js';

// The selection bar's actions follow the viewer's features.

describe('the selection bar', () => {
  it('offers Component only where the viewer has the component area', async () => {
    const v = await headlessViewer();
    const person = v.api.as({ kind: 'person' });
    const ids = (api: typeof person) => selectionModel(v.sources.interaction, api).actions.value.map((a) => a.id);
    expect(ids(person)).toContain('component');
    const without = new Proxy(person, {
      get: (target, key) =>
        key === 'area' ? (id: string) => (id === 'component' ? undefined : target.area(id)) : (target as never)[key],
    });
    expect(ids(without)).not.toContain('component');
    expect(ids(without)).toContain('focus');
  });

  it("lists Focus, Show in structure tree, Make a component, the extensions' rows by order, and the panel last", async () => {
    const v = await headlessViewer();
    const person = v.api.as({ kind: 'person' });
    const extra = observableState<readonly SelectionAction[]>([]);
    const model = selectionModel(v.sources.interaction, person, {
      openPanel: () => {},
      showInTree: () => {},
      extra: extra.asReadonly,
    });
    expect(model.actions.value.map((a) => a.id)).toEqual(['focus', 'show-in-tree', 'component', 'more']);
    expect(model.actions.value.every((a) => !!a.hint)).toBe(true);
    extra.set([{ id: 'clip', label: 'Clip…', icon: 'crop', order: 10, run: () => {} }]);
    expect(model.actions.value.map((a) => a.id)).toEqual(['focus', 'show-in-tree', 'component', 'clip', 'more']);
  });
});
