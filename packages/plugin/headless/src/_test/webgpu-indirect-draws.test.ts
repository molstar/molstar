/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { afterAll, describe, expect, it } from 'vitest';
import { Renderer } from '@molstar/graphics/gl/renderer';
import { createFrame } from '@molstar/graphics/gl/renderable';
import { IndirectRecordBytes, type MultiDrawBaseData } from '@molstar/graphics/gl/render-item';
import { type GraphicsRenderObject, isMergedRenderObject } from '@molstar/graphics/gl/render-object';
import { Scene } from '@molstar/graphics/gl/scene';
import { Camera } from '@molstar/graphics/canvas3d/camera';
import { Viewport } from '@molstar/graphics/canvas3d/camera/util';
import { HeadlessPluginContext } from '@molstar/plugin-headless/context';
import { DefaultPluginSpec } from '@molstar/plugin/default-spec';
import { createWebGPUPlugin } from './webgpu-plugin.js';
import { WebGPU } from './webgpu-module.js';

/**
 * `chains` chains of `perChain` atoms side by side: a unit each, merged; with `copies`, an assembly of that many copies
 * of each chain, a chain's at translations of its own (its segment's instances: drawn from another segment's first
 * instance, a chain shows at the other's places)
 */
function chainsPdb(chains: number, perChain: number, copies = 0) {
  const f = (v: number) => v.toFixed(3).padStart(8);
  const atoms = Array.from({ length: chains * perChain }, (_, i) => {
    const serial = String(i + 1).padStart(5);
    const c = Math.floor(i / perChain);
    const k = i % perChain;
    const residue = String(k + 1).padStart(4);
    const [x, y, z] = [(k % 3) * 3.8 + c * 14, Math.floor(k / 3) * 3.8, (i % 2) * 1.5];
    return `ATOM  ${serial}  CA  ALA ${'ABCD'[c]}${residue}    ${f(x)}${f(y)}${f(z)}  1.00  0.00           C`;
  });
  const biomt: string[] = copies > 0 ? ['REMARK 350 BIOMOLECULE: 1'] : [];
  for (let c = 0; c < chains && copies > 0; ++c) {
    biomt.push(`REMARK 350 APPLY THE FOLLOWING TO CHAINS: ${'ABCD'[c]}`);
    for (let k = 0; k < copies; ++k) {
      const n = String(c * copies + k + 1).padStart(4);
      biomt.push(`REMARK 350   BIOMT1 ${n}  1.000000  0.000000  0.000000        0.00000`);
      biomt.push(`REMARK 350   BIOMT2 ${n}  0.000000  1.000000  0.000000 ${f((c * copies + k) * 20)}`);
      biomt.push(`REMARK 350   BIOMT3 ${n}  0.000000  0.000000  1.000000        0.00000`);
    }
  }
  return `${[...biomt, ...atoms].join('\n')}\nEND\n`;
}

const Size = 64;

/** one plugin for the file's tests: its device's pipelines compiled once, each test's structure cleared after it */
let shared: Promise<HeadlessPluginContext> | undefined;
const sharedPlugin = () =>
  (shared ??= createWebGPUPlugin(DefaultPluginSpec(), { width: Size, height: Size }).then(async (plugin) => {
    await plugin.init();
    return plugin;
  }));
afterAll(async () => (await shared)?.dispose());

/** `list` with arrays of its own (the renderable may write its lists again) */
function copyList(list: MultiDrawBaseData): MultiDrawBaseData {
  return {
    ...list,
    firsts: list.firsts.slice(),
    counts: list.counts.slice(),
    offsets: list.offsets.slice(),
    instanceCounts: list.instanceCounts.slice(),
    baseVertices: list.baseVertices.slice(),
    baseInstances: list.baseInstances.slice(),
  };
}

type Representation = Parameters<
  HeadlessPluginContext['builders']['structure']['representation']['addRepresentation']
>[1];

/**
 * The opaque pass of the merged objects of `reprs` (of `chains` chains) at `size`, drawn by their own lists, the
 * segments' direct lists and the same as indirect arguments
 */
async function drawLists(
  plugin: HeadlessPluginContext,
  options: { pdb: string; chains: number; reprs: Representation[]; size: number; assembly?: boolean },
) {
  const { pdb, chains, reprs, size, assembly } = options;
  const graphics = plugin.canvas3d!.graphics;
  const data = await plugin.builders.data.rawData({ data: pdb, label: 'test' });
  const trajectory = await plugin.builders.structure.parseTrajectory(data, 'pdb');
  const model = await plugin.builders.structure.createModel(trajectory);
  const structure = await plugin.builders.structure.createStructure(
    model,
    assembly ? { name: 'assembly', params: { id: '1' } } : undefined,
  );
  for (const r of reprs) await plugin.builders.structure.representation.addRepresentation(structure, r);
  plugin.canvas3d!.commit(true);

  const objects = plugin.canvas3d!.getRenderObjects().filter((o) => isMergedRenderObject(o));
  const scene = Scene.create(graphics);
  for (const o of objects) scene.add(o);
  scene.commit();
  const camera = new Camera();
  Viewport.set(camera.viewport, 0, 0, size, size);
  const { center, radius } = scene.boundingSphere;
  camera.setState(camera.getFocus(center, radius), 0);
  camera.update();

  // drawn by the lists of each frame: recorded draw lists (render bundles) would replay the first frames' own
  const renderer = Renderer.create(graphics, { drawLists: false });
  const target = graphics.createRenderTarget(size, size, 'depth-stencil');
  await graphics.resources.preparePipelines(
    renderer.collect({ name: 'opaque' }, scene, graphics.describeTarget(target)),
  );

  // per render item the lists it draws with, taken from a frame of its own (as the renderable builds them), then drawn
  // again: as they are, or as indirect arguments of the same draws
  const captured = new Map<number, MultiDrawBaseData[]>();
  let lists: Map<number, MultiDrawBaseData[]> | undefined;
  for (const r of scene.renderables) {
    const item = r.renderItem;
    const render = item.render;
    item.render = (variant, globals, list) => {
      const own = lists?.get(item.id);
      if (!lists && list) captured.set(item.id, list.map(copyList));
      render(variant, globals, own ?? list);
    };
  }

  const frame = async (using?: Map<number, MultiDrawBaseData[]>) => {
    lists = using;
    graphics.beginRenderPass(target);
    renderer.clear(true);
    renderer.update(camera, scene, createFrame());
    renderer.setViewport(0, 0, size, size);
    renderer.renderOpaque(scene, camera);
    graphics.flush();
    return Array.from(await graphics.readRenderTarget(target, 0, 0, size, size));
  };
  // the first frame differs (what it draws with settles in it): from the second on, frames are the same
  await frame();
  const own = await frame();

  const indirect = new Map<number, MultiDrawBaseData[]>();
  const compute = graphics.compute!;
  const buffers: { destroy(): void }[] = [];
  for (const r of scene.renderables) {
    const levels = captured.get(r.renderItem.id);
    expect(levels?.length).toBeGreaterThan(0);
    // a draw a segment (a chain), its instances from the segment's first
    expect(levels![0].count).toBe(chains);
    const o = objects.find((o) => o.id === r.id)!;
    const indexed = o.type === 'mesh';
    const stride = indexed ? IndirectRecordBytes.indexed : IndirectRecordBytes.vertex;
    indirect.set(
      r.renderItem.id,
      levels!.map((list) => {
        const args = new Uint32Array((list.count * stride) / 4);
        for (let i = 0; i < list.count; ++i) {
          const record = indexed
            ? [
                list.counts[i],
                list.instanceCounts[i],
                list.offsets[i] / 4,
                list.baseVertices[i] ?? 0,
                list.baseInstances[i],
              ]
            : [list.counts[i], list.instanceCounts[i], list.firsts[i], list.baseInstances[i]];
          args.set(record, (i * stride) / 4);
        }
        // at an offset: the records start where the list says
        const byteOffset = 64;
        const buffer = compute.createBuffer(byteOffset + Math.max(4, args.byteLength), 'test indirect args');
        buffer.write(args, byteOffset);
        buffers.push(buffer);
        return { ...list, indirect: { buffer, byteOffset } };
      }),
    );
  }
  const directImage = await frame(captured);
  const indirectImage = await frame(indirect);
  lists = undefined;
  for (const b of buffers) b.destroy();
  renderer.dispose();
  scene.clear();
  return { types: objects.map((o) => o.type).sort(), own, directImage, indirectImage };
}

describe.skipIf(!WebGPU.module)('indirect draw lists on WebGPU', () => {
  // quick: sphere impostors (`drawIndirect`) of two chains of three atoms, 16 px
  it('draw what the same direct lists draw, from instances other than the first', async () => {
    const plugin = await sharedPlugin();
    try {
      expect(plugin.canvas3d!.graphics.capabilities.indirectInstanceRanges).toBe(true);
      const { types, own, directImage, indirectImage } = await drawLists(plugin, {
        // two copies a chain, each chain's at places of its own: a segment drawn from another one's first instance
        // shows (its spheres at the other's places)
        pdb: chainsPdb(2, 3, 2),
        chains: 2,
        reprs: [{ type: 'spacefill', typeParams: {}, color: 'chain-id' }],
        size: 16,
        assembly: true,
      });
      expect(types).toEqual(['spheres']);
      expect(new Set(own).size).toBeGreaterThan(2);
      expect(directImage).toEqual(own);
      expect(indirectImage).toEqual(own);
    } finally {
      await plugin.clear();
    }
  });

  // Medium: the pipelines of two representations compiled on a new device, a second or more
  it(
    'draw what the same direct lists draw, indexed or not, from instances other than the first',
    { tags: ['medium'] },
    async () => {
      const plugin = await sharedPlugin();
      try {
        expect(plugin.canvas3d!.graphics.capabilities.indirectInstanceRanges).toBe(true);
        const { types, own, directImage, indirectImage } = await drawLists(plugin, {
          // two copies a chain, each chain's at places of its own: a segment drawn from another one's first instance
          // shows (its spheres and meshes at the other's places)
          pdb: chainsPdb(4, 12, 2),
          chains: 4,
          assembly: true,
          // sphere impostors (`drawIndirect`) and meshes (`drawIndexedIndirect`)
          reprs: [
            { type: 'spacefill', typeParams: {}, color: 'chain-id' },
            { type: 'ball-and-stick', typeParams: { sizeFactor: 0.2, tryUseImpostor: false }, color: 'chain-id' },
          ],
          size: Size,
        });
        expect(types).toEqual(['mesh', 'spheres']);
        // the objects are drawn (not only the background), the lists draw what the objects' own draw
        expect(new Set(own).size).toBeGreaterThan(2);
        expect(directImage).toEqual(own);
        expect(indirectImage).toEqual(own);
      } finally {
        await plugin.clear();
      }
    },
    60_000,
  );
});
