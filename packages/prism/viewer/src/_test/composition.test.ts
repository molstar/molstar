/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { afterEach, describe, expect, it } from 'vitest';
import { ApiResult } from '../api/area.js';
import { DefaultPrismViewerSpec } from '../default-spec.js';
import { AnnotationsFeature } from '../features/annotations/index.js';
import { HomeFeature } from '../features/home/index.js';
import { StructureFeature } from '../features/structure/index.js';
import { StructurePanelFeature } from '../features/structure-panel/index.js';
import { createPrismViewerHeadless } from '../headless.js';
import type { ViewerExtension } from '../spec.js';
import type { PrismViewer } from '../viewer.js';

// Composing a viewer (viewer-spec.md §3.1): the built-in features' areas are there when the spec lists the feature;
// left out, reading one says which feature to add, its operations are not defined, and what needs it is refused or
// goes without it.

const viewers: PrismViewer[] = [];
afterEach(async () => {
  for (const v of viewers.splice(0)) await v.dispose();
});

async function viewer(extensions: ViewerExtension[]) {
  const v = await createPrismViewerHeadless({ spec: { ...DefaultPrismViewerSpec(), extensions } });
  viewers.push(v);
  return v;
}

describe('a viewer composed of fewer features', () => {
  it('has the core areas; a built-in area left out names its feature, and its operations are not defined', async () => {
    const v = await viewer([]);
    expect(v.api.camera).toBeDefined();
    expect('volume' in v.api).toBe(false);
    expect(v.api.area('volume')).toBeUndefined();
    expect(() => v.api.volume).toThrow(
      "The 'volume' area is not in this viewer: add VolumeFeature to its spec's extensions (DefaultPrismViewerSpec has it).",
    );
    expect(() => v.api.representation).toThrow(/add StructureFeature/);
    expect(v.api.describe('structure.preset')).toBeUndefined();
    expect(await v.api.run('structure.preset', {})).toMatchObject({ ok: false, code: 'not-found' });
  });

  it('has what a listed feature brings', async () => {
    const v = await viewer([StructureFeature]);
    expect(v.api.area('structure')).toBeDefined();
    expect(v.api.describe('representation.color')).toBeDefined();
    expect(v.managers.commands.get('structure.preset')).toBeDefined();
    expect(v.api.area('trajectory')).toBeUndefined();
  });

  it('refuses a feature listed without the ones it needs before it', async () => {
    await expect(viewer([AnnotationsFeature])).rejects.toThrow(
      "extension 'annotations' needs 'structure' before it in the spec's extensions",
    );
    await expect(viewer([AnnotationsFeature, StructureFeature])).rejects.toThrow(/needs 'structure'/);
  });

  it("activates an extension's annotations only where the annotations feature is", async () => {
    const scores: ViewerExtension = {
      id: 'scores',
      title: 'Scores',
      annotations: [
        {
          id: 'test-score',
          title: 'Test score',
          isApplicable: () => true,
          apply: async () => ApiResult.ok('Shown.', {}),
        },
      ],
    };
    const without = await viewer([StructureFeature, scores]);
    expect(without.api.area('annotations')).toBeUndefined();
    const withIt = await viewer([StructureFeature, AnnotationsFeature, scores]);
    expect(withIt.api.annotations.providers.value.map((p) => p.id)).toContain('test-score');
  });

  it("makes a panel feature's view models with it; without it, reading one names the feature", async () => {
    const bare = await viewer([StructureFeature]);
    expect('scene' in bare.views).toBe(false);
    expect(() => bare.views.scene).toThrow(
      "The 'scene' view model is not in this viewer: add StructurePanelFeature to its spec's extensions.",
    );
    expect(() => bare.views.homeActions).toThrow(/add HomeFeature/);
    const full = await viewer([StructureFeature, StructurePanelFeature, HomeFeature]);
    expect(full.views.scene.roots.value).toEqual([]);
    expect(full.views.homeActions.value.length).toBeGreaterThan(0);
    await expect(viewer([StructurePanelFeature])).rejects.toThrow(/'structure-panel' needs 'structure'/);
  });
});
