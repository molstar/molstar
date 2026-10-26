/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { describe, expect, it } from 'vitest';
import { ValueCell } from '@molstar/core/util';
import {
  createSignatureProgramItem,
  getObjectSignature,
  getSignatureKey,
  type ObjectSignature,
} from '@molstar/graphics/gl/precompile';
import { Renderer } from '@molstar/graphics/gl/renderer';
import { Scene } from '@molstar/graphics/gl/scene';
import type { GlobalDefines } from '@molstar/graphics/gl/renderable/schema';
import { HeadlessPluginContext } from '@molstar/plugin-headless/context';
import { DefaultPluginSpec } from '@molstar/plugin/default-spec';
import { WebGPU } from './webgpu-module.js';

const Pdb = `ATOM      1  N   ALA A   1      -1.458   0.000   0.000  1.00  0.00           N
ATOM      2  CA  ALA A   1       0.000   0.000   0.000  1.00  0.00           C
ATOM      3  C   ALA A   1       0.551   1.420   0.000  1.00  0.00           C
ATOM      4  O   ALA A   1      -0.200   2.390   0.000  1.00  0.00           O
ATOM      5  CB  ALA A   1       0.540  -0.760  -1.210  1.00  0.00           C
ATOM      6  N   GLY A   2       1.870   1.550   0.000  1.00  0.00           N
ATOM      7  CA  GLY A   2       2.520   2.860   0.000  1.00  0.00           C
ATOM      8  C   GLY A   2       4.030   2.720   0.000  1.00  0.00           C
ATOM      9  O   GLY A   2       4.570   1.610   0.000  1.00  0.00           O
HETATM   10  O   HOH A   3       6.000   0.000   0.000  1.00  0.00           O
END
`;

const Passes: Renderer.Pass[] = [
  { name: 'opaque' },
  { name: 'wboitTransparent' },
  { name: 'blendedTransparent' },
  { name: 'depthOpaque' },
  { name: 'depthTransparent' },
  { name: 'markingDepth' },
  { name: 'markingMask' },
  { name: 'pick', variant: 'pick' },
];

const createPlugin = async () => {
  const plugin = await HeadlessPluginContext.create(
    { webgpu: WebGPU.module! },
    DefaultPluginSpec(),
    { width: 32, height: 32 },
    { backend: 'webgpu' },
  );
  await plugin.init();
  const graphics = plugin.canvas3d!.graphics;
  // the draw pass's formats, with stencil
  const target = graphics.describeTarget(graphics.createRenderTarget(8, 8, 'depth-stencil'));
  return { plugin, graphics, renderer: Renderer.create(graphics), target };
};

/** meshes, sphere and cylinder impostors, text; opaque, transparent, double sided */
async function addRepresentations(plugin: HeadlessPluginContext) {
  const data = await plugin.builders.data.rawData({ data: Pdb, label: 'test' });
  const trajectory = await plugin.builders.structure.parseTrajectory(data, 'pdb');
  const model = await plugin.builders.structure.createModel(trajectory);
  const structure = await plugin.builders.structure.createStructure(model);
  const { representation } = plugin.builders.structure;
  await representation.addRepresentation(structure, { type: 'cartoon' });
  await representation.addRepresentation(structure, { type: 'ball-and-stick', typeParams: { alpha: 0.5 } });
  await representation.addRepresentation(structure, { type: 'label' });
  await representation.addRepresentation(structure, {
    type: 'molecular-surface',
    typeParams: { doubleSided: true, transparentBackfaces: 'on' },
  });
  plugin.canvas3d!.commit(true);
}

async function addScene(plugin: HeadlessPluginContext) {
  await addRepresentations(plugin);
  // a scene of the plugin's render objects, as Canvas3D draws them
  const scene = Scene.create(plugin.canvas3d!.graphics, Transparency, Globals);
  for (const ro of plugin.canvas3d!.getRenderObjects()) scene.add(ro);
  scene.commit();
  return scene;
}

const Transparency = 'wboit';
const Globals: GlobalDefines = { dColorMarker: true, dLightCount: 1 };

describe.skipIf(!WebGPU.module)('object signatures on WebGPU', () => {
  it('prepare the programs and pipelines of their objects, also in another context', async () => {
    const a = await createPlugin();
    const b = await createPlugin();
    try {
      const scene = await addScene(a.plugin);
      expect(scene.renderables.length).toBeGreaterThan(3);
      const { resources } = a.graphics;
      for (const pass of Passes) {
        await resources.preparePipelines(a.renderer.collect(pass, scene, a.target, { all: true }));
      }

      // signatures as plain data (as a precompile list stores them), deduplicated
      const signatures = new Map<string, ObjectSignature>();
      for (const r of scene.renderables) {
        const signature = getObjectSignature(r);
        expect(signature).toBeDefined();
        signatures.set(getSignatureKey(signature!), JSON.parse(JSON.stringify(signature)));
      }
      const items = [...signatures.values()].map((signature) => ({
        signature,
        item: createSignatureProgramItem(a.graphics, signature, Transparency, Globals)!,
      }));

      // in the same context: the objects' programs and pipelines, nothing new
      const programs = a.graphics.stats.resourceCounts.program;
      const tracked = resources.pending.trackedCount;
      let requestCount = 0;
      for (const pass of Passes) {
        const requests = a.renderer.collectSignatures(pass, items, a.target, { all: true });
        expect(requests.length).toBeGreaterThan(0);
        requestCount += requests.length;
        await resources.preparePipelines(requests);
      }
      expect(requestCount).toBeGreaterThanOrEqual(items.length * Passes.length);
      expect(resources.pending.trackedSince(tracked)).toEqual([]);
      expect(a.graphics.stats.resourceCounts.program).toBe(programs);
      for (const { item } of items) item.destroy();

      // in a new context, prepared from the signatures before the objects exist: drawing them compiles nothing more
      // into the target (Canvas3D's own targets aside)
      const bItems = [...signatures.values()].map((signature) => ({
        signature,
        item: createSignatureProgramItem(b.graphics, signature, Transparency, Globals)!,
      }));
      const bPending = b.graphics.resources.pending;
      const beforePrepare = bPending.trackedCount;
      for (const pass of Passes) {
        await b.graphics.resources.preparePipelines(
          b.renderer.collectSignatures(pass, bItems, b.target, { all: true }),
        );
      }
      expect(bPending.trackedCount).toBeGreaterThan(beforePrepare);
      const bScene = await addScene(b.plugin);
      // what Canvas3D compiles for its own scene and passes settled first
      await bPending.whenReady();
      const before = bPending.trackedCount;
      for (const pass of Passes) {
        await b.graphics.resources.preparePipelines(b.renderer.collect(pass, bScene, b.target, { all: true }));
      }
      // pipeline labels end with their key, `<render state>|<target>|<draw mode>`
      const intoTarget = bPending
        .trackedSince(before)
        .filter((label) => label.split(' ').at(-1)!.split('|').slice(1, -1).join('|') === b.target.key);
      expect(intoTarget).toEqual([]);

      // control: without the signatures prepared, the same compiles something into the target
      const c = await createPlugin();
      try {
        const cPending = c.graphics.resources.pending;
        const cScene = await addScene(c.plugin);
        await cPending.whenReady();
        const cBefore = cPending.trackedCount;
        for (const pass of Passes) {
          await c.graphics.resources.preparePipelines(c.renderer.collect(pass, cScene, c.target, { all: true }));
        }
        const cInto = cPending
          .trackedSince(cBefore)
          .filter((label) => label.split(' ').at(-1)!.split('|').slice(1, -1).join('|') === c.target.key);
        expect(cInto.length).toBeGreaterThan(0);
      } finally {
        c.plugin.dispose();
      }
    } finally {
      a.plugin.dispose();
      b.plugin.dispose();
    }
  });

  it('an updated representation records its objects again, not the rest of the scene', async () => {
    const { plugin } = await createPlugin();
    try {
      const data = await plugin.builders.data.rawData({ data: Pdb, label: 'test' });
      const trajectory = await plugin.builders.structure.parseTrajectory(data, 'pdb');
      const model = await plugin.builders.structure.createModel(trajectory);
      const structure = await plugin.builders.structure.createStructure(model);
      const { representation } = plugin.builders.structure;
      const cartoon = (await representation.addRepresentation(structure, { type: 'cartoon' })).data!.repr;
      const spacefill = (await representation.addRepresentation(structure, { type: 'spacefill' })).data!.repr;
      const canvas = plugin.canvas3d!;
      canvas.commit(true);

      // the types of the signatures recorded since `keys`
      const keysOf = () => new Map(canvas.getPrecompileList().objects.map((o) => [getSignatureKey(o), o.type]));
      const newTypes = (keys: Map<string, string>) =>
        new Set([...keysOf()].filter(([k]) => !keys.has(k)).map(([, type]) => type));
      const before = keysOf();

      // a signature value changed on the objects of both (meshes, spheres), an update of the cartoon only
      for (const repr of [cartoon, spacefill]) {
        for (const ro of repr.renderObjects)
          ValueCell.update(ro.values.uDoubleSided!, !ro.values.uDoubleSided!.ref.value);
      }
      canvas.update(cartoon);
      expect(newTypes(before)).toEqual(new Set(['mesh']));

      // the spheres' change is recorded with their own update
      canvas.update(spacefill);
      expect(newTypes(before)).toEqual(new Set(['mesh', 'spheres']));
    } finally {
      plugin.dispose();
    }
  });

  it('a recorded precompile list prepares what frames draw the objects with, before they exist', async () => {
    const a = await createPlugin();
    const b = await createPlugin();
    const c = await createPlugin();
    try {
      const canvasA = a.plugin.canvas3d!;
      canvasA.resetPrecompileList();
      await addRepresentations(a.plugin);
      const list = JSON.parse(JSON.stringify(canvasA.getPrecompileList()));
      expect(list.objects.length).toBeGreaterThan(3);

      // the context's own (empty scene) settled, then the list for the current transparency (wboit) and blended
      const canvasB = b.plugin.canvas3d!;
      const pending = b.graphics.resources.pending;
      await canvasB.prepare({ interaction: true });
      await pending.whenReady();
      const beforePrecompile = pending.trackedCount;
      await canvasB.precompile(list, [{}, { transparency: 'blended' }]);
      expect(pending.trackedCount).toBeGreaterThan(beforePrecompile);

      // the objects added (only to Canvas3D, not to a scene of the test): their frames and the first interaction compile
      // nothing more
      const before = pending.trackedCount;
      await addRepresentations(b.plugin);
      await canvasB.prepare({ interaction: true });
      await pending.whenReady();
      expect(pending.trackedSince(before)).toEqual([]);

      // control: without the list, the same compiles
      const canvasC = c.plugin.canvas3d!;
      const cPending = c.graphics.resources.pending;
      await canvasC.prepare({ interaction: true });
      await cPending.whenReady();
      const cBefore = cPending.trackedCount;
      await addRepresentations(c.plugin);
      await canvasC.prepare({ interaction: true });
      await cPending.whenReady();
      expect(cPending.trackedSince(cBefore).length).toBeGreaterThan(0);
    } finally {
      a.plugin.dispose();
      b.plugin.dispose();
      c.plugin.dispose();
    }
  });
});
