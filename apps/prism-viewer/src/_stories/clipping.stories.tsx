/** @jsxImportSource @molstar/prism-ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { Meta, StoryObj } from '@storybook/web-components-vite';
import { ClippingExtension, ClippingSectionId } from '@molstar/clipping-extension/prism-viewer';
import { story } from '@molstar/prism-ui/testing/index';
import type { PrismViewer } from '@molstar/prism-viewer';
import type { PrismViewerSpec } from '@molstar/prism-viewer/spec';
import { ExampleFiles, RealViewerStory } from '@molstar/prism-viewer/testing/real-viewer';

// The clipping extension in the Prism viewer: cuts through representations (a plane, a sphere, a box, a cylinder, a
// cone), made and edited in the Selection tab's Clipping section, from the selection menu's Clip…, or by the API.
// 1CBS (data/examples): its polymer drawn as a molecular surface over the cartoon, retinoic acid inside. No network.

const meta: Meta = { title: 'Viewer/With extensions/Clipping', parameters: { layout: 'fullscreen' } };
export default meta;

const withClipping = (base: PrismViewerSpec): PrismViewerSpec => ({
  ...base,
  extensions: [...(base.extensions ?? []), ClippingExtension],
});

const ok = (r: { ok: boolean; message?: string }) => {
  if (!r.ok) throw new Error(r.message);
};

/** 1CBS with its polymer as a surface too, the polymer selected, the Clipping section open; then the story's steps */
function clippingStory(steps?: (v: PrismViewer) => Promise<void>): StoryObj {
  return {
    render: () =>
      story(() => (
        <RealViewerStory
          spec={withClipping}
          load={[ExampleFiles.cbs]}
          setup={async (v) => {
            const polymer = v.plugin.managers.structure.hierarchy.current.structures[0].components.find(
              (c) => c.cell.obj?.label === 'Polymer',
            );
            if (polymer)
              ok(await v.api.representation.set({ component: polymer.cell.transform.ref, type: 'molecular-surface' }));
            ok(await v.api.select.set({ target: 'polymer' }));
            await v.api.view.open({ tab: 'selection', section: ClippingSectionId });
            await steps?.(v);
          }}
        />
      )),
  };
}

/** The section, nothing cut yet: a cut by shape, narrowed to the types the selection is drawn as */
export const Section: StoryObj = clippingStory();

/** A plane through the surface only: the cartoon and the ligand inside show where the surface is cut away */
export const PlaneThroughTheSurface: StoryObj = clippingStory(async (v) => {
  ok(await v.api.clipping!.add({ shape: 'plane', types: ['molecular-surface'] }));
});

/** A sphere window into the surface, and its numbers open */
export const SphereWindow: StoryObj = clippingStory(async (v) => {
  ok(await v.api.clipping!.add({ shape: 'sphere', types: ['molecular-surface'] }));
});

/** Every shape, one cut each on the surface */
export const EveryShape: StoryObj = clippingStory(async (v) => {
  for (const shape of ['plane', 'sphere', 'cube', 'cylinder', 'infiniteCone'] as const) {
    ok(await v.api.clipping!.add({ shape, types: ['molecular-surface'] }));
  }
});
