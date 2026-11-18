/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The bounds as the viewer's composition holds them (_spec/README.md "Configuration"): items with the security
// table's defaults, the patch's nesting check and the store's table bounds reading them.

import type { Patch } from '@molstar/core/repr/patch';
import { describe, expect, it } from 'vitest';
import { DefaultViewerConfig, viewerConfig, viewerConfigOf, viewerItem } from '../config.js';
import { SpecModel } from '../models/spec.js';
import { StoreModel } from '../models/store.js';
import { proteinOnGene } from './fixture.js';

const hideLetters = (): Patch => ({
  ops: [{ op: 'set', at: ['views', 'main', 'protein', 'letters', 'l'], key: 'hidden', value: true }],
});

describe('the bounds items', () => {
  it('have the security table’s defaults', () => {
    expect(DefaultViewerConfig.bounds).toEqual({
      maxBytes: 64 * 1024 * 1024,
      maxDepth: 256,
      maxQueryDepth: 64,
      maxColumnLength: 2 ** 31 - 1,
      maxNodeIssues: 100,
      maxPlaces: 10,
      maxBlockValues: 2 ** 24,
    });
  });

  it('are set as items of a viewer spec', () => {
    const config = viewerConfigOf([
      [viewerItem('bounds.maxDepth'), 4],
      [viewerItem('bounds.maxBlockValues'), 100],
    ]);
    expect(config.bounds).toMatchObject({ maxDepth: 4, maxBlockValues: 100, maxBytes: 64 * 1024 * 1024 });
  });
});

describe('a patch under the nesting bound', () => {
  it('is applied at the bound and refused one level past it, as the whole spec’s check does', () => {
    // spec / views / main / protein / letters / l: the layer is the sixth level, the root the zeroth
    const at = new SpecModel(
      proteinOnGene(),
      {},
      viewerConfig({ bounds: { ...DefaultViewerConfig.bounds, maxDepth: 5 } }),
    );
    const ok = at.applyPatch(hideLetters());
    expect(ok.ok && ok.changed).toBe(true);
    const past = new SpecModel(
      proteinOnGene(),
      {},
      viewerConfig({ bounds: { ...DefaultViewerConfig.bounds, maxDepth: 4 } }),
    );
    const before = past.head.value.revision;
    const refused = past.applyPatch(hideLetters());
    expect(refused.ok).toBe(false);
    expect(refused.issues.map((i) => i.code)).toContain('seq.lang.bound.depth');
    expect(past.head.value.revision).toBe(before);
  });

  it('is applied by default', () => {
    const spec = new SpecModel(proteinOnGene());
    expect(spec.applyPatch(hideLetters()).ok).toBe(true);
  });
});

describe('the store’s table bounds', () => {
  it('are the composition’s', () => {
    const store = new StoreModel(viewerConfig({ bounds: { ...DefaultViewerConfig.bounds, maxBlockValues: 7 } }));
    expect(store.store.settings).toMatchObject({ maxBlockValues: 7, maxColumnLength: 2 ** 31 - 1 });
  });
});
