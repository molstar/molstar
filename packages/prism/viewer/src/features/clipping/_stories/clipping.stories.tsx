/** @jsxImportSource @molstar/prism-ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Ludovic Autin <ludovic.autin@gmail.com>
 */

import type { Meta, StoryObj } from '@storybook/web-components-vite';
import { ClippingSectionId } from '../index.js';
import { story } from '@molstar/prism-ui/testing/index';
import type { PrismViewer } from '../../../index.js';
import { ExampleFiles, RealViewerStory } from '../../../testing/real-viewer.js';

// The clipping feature (clipping.md): cuts through representations (a plane, a sphere, a box, a cylinder, a
// cone), made and edited in the Selection tab's Clipping section, from the selection menu's Clip…, or by the API.
// 1CBS (data/examples): its polymer drawn as a molecular surface over the cartoon, retinoic acid inside. No network.

const meta: Meta = { title: 'Viewer/Clipping', parameters: { layout: 'fullscreen' } };
export default meta;

const ok = (r: { ok: boolean; message?: string }) => {
  if (!r.ok) throw new Error(r.message);
};

/** 1CBS with its polymer as a surface too, the polymer selected, the Clipping section open; then the story's steps */
function clippingStory(steps?: (v: PrismViewer) => Promise<void>): StoryObj {
  return {
    render: () =>
      story(() => (
        <RealViewerStory
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
  ok(await v.api.clipping.add({ shape: 'plane', types: ['molecular-surface'] }));
});

/** A sphere window into the surface, and its numbers open */
export const SphereWindow: StoryObj = clippingStory(async (v) => {
  ok(await v.api.clipping.add({ shape: 'sphere', types: ['molecular-surface'] }));
});

/** Every shape, one cut each on the surface */
export const EveryShape: StoryObj = clippingStory(async (v) => {
  for (const shape of ['plane', 'sphere', 'cube', 'cylinder', 'infiniteCone'] as const) {
    ok(await v.api.clipping.add({ shape, types: ['molecular-surface'] }));
  }
});

/** A plane with its handles in 3D: drag the orange knob to move it, the blue arrow to slide it, the green ring to turn it */
export const PlaneHandles: StoryObj = clippingStory(async (v) => {
  ok(await v.api.clipping.add({ shape: 'plane', types: ['molecular-surface'] }));
  v.api.clipping.gizmos.handleOnly(1);
});

/** A box with its handles: the purple corner knob sizes it */
export const BoxHandles: StoryObj = clippingStory(async (v) => {
  ok(await v.api.clipping.add({ shape: 'cube', types: ['molecular-surface'] }));
  v.api.clipping.gizmos.handleOnly(1);
});

/** Clip this… from the structure tree's row menu, on the polymer's surface: a plane through it, active */
export const ClipThisFromTheTree: StoryObj = clippingStory(async (v) => {
  const polymer = v.plugin.managers.structure.hierarchy.current.structures[0].components.find(
    (c) => c.cell.obj?.label === 'Polymer',
  );
  const surface = polymer?.representations.find(
    (r) => (r.cell.transform.params as { type?: { name: string } }).type?.name === 'molecular-surface',
  );
  if (!surface) throw new Error('No surface');
  const done = await v.api.clipping.clipThis(v.api.as({ kind: 'person' }), surface.cell.transform.ref);
  if (done && done.cut) v.api.clipping.gizmos.handleOnly(done.cut.index);
});

/** A sphere around the ligand inside the surface (centered on what it targets, with room) */
export const AroundTheLigand: StoryObj = clippingStory(async (v) => {
  ok(await v.api.select.set({ target: 'ligand' }));
  const polymer = v.plugin.managers.structure.hierarchy.current.structures[0].components.find(
    (c) => c.cell.obj?.label === 'Polymer',
  );
  const surface = polymer?.representations.find(
    (r) => (r.cell.transform.params as { type?: { name: string } }).type?.name === 'molecular-surface',
  );
  if (!surface) throw new Error('No surface');
  ok(await v.api.clipping.add({ shape: 'sphere', node: `#${surface.cell.transform.ref}`, target: 'ligand' }));
  v.api.clipping.gizmos.handleOnly(1);
});

/** Every shape drawn in the scene (clip volumes, saved with it): translucent walls, solid edges, no handles */
export const DrawnInTheScene: StoryObj = clippingStory(async (v) => {
  for (const shape of ['plane', 'sphere', 'cube', 'cylinder', 'infiniteCone'] as const) {
    ok(await v.api.clipping.add({ shape, types: ['molecular-surface'] }));
  }
  for (let i = 1; i <= 5; i++) ok(await v.api.clipping.show({ cut: i, on: true }));
});
