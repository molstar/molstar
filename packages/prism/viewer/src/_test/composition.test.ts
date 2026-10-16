/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { afterEach, describe, expect, it, vi } from 'vitest';
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

  it('names the feature of a built-in operation in a core area, or a built-in area, it left out', async () => {
    const v = await viewer([StructureFeature]);
    expect(v.api.describe('load.entry')).toBeUndefined();
    expect(() => v.api.load.entry).toThrow(
      "The operation 'load.entry' is not in this viewer: add LoadersFeature to its spec's extensions (DefaultPrismViewerSpec has it).",
    );
    expect(() => v.api.canvas.lighting).toThrow(/add CanvasSettingsFeature/);
    expect(() => v.api.session.save).toThrow(/add ExportSessionFeature/);
    expect(() => v.api.snapshot).toThrow(/add SnapshotsFeature/);
    // The core's own operations are there.
    expect(v.api.load.files).toBeTypeOf('function');
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

  it('checks overrides and keybindings against the commands once they are in; headless, the page-only ones pass', async () => {
    // Preset's and Help's commands are the page's: overriding them headless is not a mistake.
    const v = await createPrismViewerHeadless({
      spec: { ...DefaultPrismViewerSpec(), overrides: { 'help.show': null, 'core.preset.auto': null } },
    });
    viewers.push(v);
    await expect(
      createPrismViewerHeadless({ spec: { ...DefaultPrismViewerSpec(), overrides: { 'no.such-command': null } } }),
    ).rejects.toThrow("override 'no.such-command': there is no command 'no.such-command'");
    const warn = vi.spyOn(console, 'warn').mockImplementation(() => {});
    const k = await createPrismViewerHeadless({
      spec: { ...DefaultPrismViewerSpec(), ui: { keybindings: [{ command: 'no.such-command', key: { key: 'x' } }] } },
    });
    viewers.push(k);
    expect(warn).toHaveBeenCalledWith("ui.keybindings: there is no command 'no.such-command' (its key is not bound).");
    warn.mockRestore();
  });

  it('says on viewer.events.error that an extension did not load', async () => {
    const v = await viewer([StructureFeature]);
    const errors: string[] = [];
    const sub = v.events.error.subscribe((e) => errors.push(`${e.source}: ${e.message}`));
    const broken = {
      id: 'broken',
      title: 'Broken',
      stub: { parts: [] },
      load: async () => {
        throw new Error('no chunk');
      },
    };
    await v.extensions.load('broken').catch(() => {});
    // Not listed: nothing to load; listed lazily, its failure is said.
    const w = await viewer([StructureFeature, broken as never]);
    const later: string[] = [];
    const sub2 = w.events.error.subscribe((e) => later.push(`${e.source}: ${e.message}`));
    await w.extensions.load('broken');
    expect(later.some((l) => l.startsWith('extension: broken could not load'))).toBe(true);
    sub.unsubscribe();
    sub2.unsubscribe();
  });
});
