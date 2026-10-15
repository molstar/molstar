/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { describe, expect, it } from 'vitest';
import { headlessViewer } from '../../../_test/harness.js';
import { DefaultPrismViewerSpec } from '../../../default-spec.js';
import { isLazyExtension, type ViewerExtension } from '../../../spec.js';
import { LazyFeatures } from '../index.js';

// The lazy features' stubs are what their features bring (viewer-spec.md §3: a load that differs fails).

describe('lazy features', () => {
  for (const lazy of Object.values(LazyFeatures)) {
    it(`${lazy.id}: the stub declares what the feature brings, and the feature loads`, async () => {
      const feature: ViewerExtension = await lazy.load();
      const v = await headlessViewer();
      expect(v.extensions.status(lazy.id)?.value).toBe('stub');
      expect(await v.extensions.load(lazy.id)).toBe(true);
      expect(v.extensions.status(lazy.id)?.value).toBe('loaded');
      for (const op of lazy.stub.operations ?? []) {
        expect(v.api.describe(op.id)?.definition.title).toBe(op.title);
      }
      for (const c of lazy.stub.commands ?? []) expect(v.managers.commands.get(c.id)?.title).toBe(c.title);
      expect(feature.builtIn).toBe(lazy.builtIn);
    });
  }

  it('are lazy in the default spec', () => {
    const ids = (DefaultPrismViewerSpec().extensions ?? []).filter(isLazyExtension).map((e) => e.id);
    expect(ids).toEqual(expect.arrayContaining(Object.keys(LazyFeatures)));
  });
});
