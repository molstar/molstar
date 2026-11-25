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
import { WebGPU } from './webgpu-module.js';

/** `chains` chains of `perChain` atoms side by side: a unit each, merged */
function chainsPdb(chains: number, perChain: number) {
  const f = (v: number) => v.toFixed(3).padStart(8);
  const atoms = Array.from({ length: chains * perChain }, (_, i) => {
    const serial = String(i + 1).padStart(5);
    const c = Math.floor(i / perChain);
    const k = i % perChain;
    const residue = String(k + 1).padStart(4);
    const [x, y, z] = [(k % 3) * 3.8 + c * 14, Math.floor(k / 3) * 3.8, (i % 2) * 1.5];
    return `ATOM  ${serial}  CA  ALA ${'ABCD'[c]}${residue}    ${f(x)}${f(y)}${f(z)}  1.00  0.00           C`;
  });
  return `${atoms.join('\n')}\nEND\n`;
}

const Size = 64;

/** one plugin for the file's tests: its device's pipelines compiled once, each test's structure cleared after it */
let shared: Promise<HeadlessPluginContext> | undefined;
const sharedPlugin = () =>
  (shared ??= HeadlessPluginContext.create(
    { webgpu: WebGPU.module! },
    DefaultPluginSpec(),
    { width: Size, height: Size },
    { backend: 'webgpu' },
  ).then(async (plugin) => {
    await plugin.init();
    return plugin;
  }));
afterAll(async () => (await shared)?.dispose());

/** A list of the segments' draws of merged `o`: a draw per segment, its instances from the segment's first */
function getSegmentList(o: GraphicsRenderObject): MultiDrawBaseData {
  if (!isMergedRenderObject(o)) throw new Error('expected a merged render object');
  const { segments } = o.merged;
  const hasElements = !!(o.values as { elements?: unknown }).elements;
  const n = segments.count;
  const list: MultiDrawBaseData = {
    firsts: new Int32Array(n),
    counts: new Int32Array(n),
    offsets: new Int32Array(n),
    instanceCounts: new Int32Array(n),
    baseVertices: new Int32Array(n),
    baseInstances: new Uint32Array(n),
    count: n,
    uniforms: [],
  };
  for (let i = 0; i < n; ++i) {
    list.firsts[i] = hasElements ? 0 : segments.drawBases[i];
    list.counts[i] = o.merged.records[i].drawCount.ref.value;
    list.offsets[i] = segments.elementOffsets[i];
    list.instanceCounts[i] = segments.instanceCounts[i];
    list.baseInstances[i] = segments.instanceBases[i];
  }
  return list;
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
  options: { pdb: string; chains: number; reprs: Representation[]; size: number },
) {
  const { pdb, chains, reprs, size } = options;
  const graphics = plugin.canvas3d!.graphics;
  const data = await plugin.builders.data.rawData({ data: pdb, label: 'test' });
  const trajectory = await plugin.builders.structure.parseTrajectory(data, 'pdb');
  const model = await plugin.builders.structure.createModel(trajectory);
  const structure = await plugin.builders.structure.createStructure(model);
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

  const renderer = Renderer.create(graphics);
  const target = graphics.createRenderTarget(size, size, 'depth-stencil');
  await graphics.resources.preparePipelines(
    renderer.collect({ name: 'opaque' }, scene, graphics.describeTarget(target)),
  );

  // per render item the lists to draw instead of its own: direct, or the same as indirect arguments
  const direct = new Map<number, MultiDrawBaseData>();
  const indirect = new Map<number, MultiDrawBaseData>();
  const compute = graphics.compute!;
  for (const r of scene.renderables) {
    const o = objects.find((o) => o.id === r.id)!;
    const list = getSegmentList(o);
    expect(list.count).toBe(chains);
    direct.set(r.renderItem.id, list);
    const indexed = o.type === 'mesh';
    const stride = indexed ? IndirectRecordBytes.indexed : IndirectRecordBytes.vertex;
    const args = new Uint32Array((list.count * stride) / 4);
    for (let i = 0; i < list.count; ++i) {
      const record = indexed
        ? [list.counts[i], list.instanceCounts[i], list.offsets[i] / 4, 0, list.baseInstances[i]]
        : [list.counts[i], list.instanceCounts[i], list.firsts[i], list.baseInstances[i]];
      args.set(record, (i * stride) / 4);
    }
    // at an offset: the records start where the list says
    const byteOffset = 64;
    const buffer = compute.createBuffer(byteOffset + args.byteLength, 'test indirect args');
    buffer.write(args, byteOffset);
    indirect.set(r.renderItem.id, { ...list, indirect: { buffer, byteOffset } });
  }
  let lists: Map<number, MultiDrawBaseData> | undefined;
  for (const r of scene.renderables) {
    const item = r.renderItem;
    const render = item.render;
    item.render = (variant, globals, list) => {
      const own = lists?.get(item.id);
      render(variant, globals, own ? [own] : list);
    };
  }

  const frame = async (using?: Map<number, MultiDrawBaseData>) => {
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
  const directImage = await frame(direct);
  const indirectImage = await frame(indirect);
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
        pdb: chainsPdb(2, 3),
        chains: 2,
        // colored by chain: a segment drawn with another one's instances (their segment data) shows
        reprs: [{ type: 'spacefill', typeParams: {}, color: 'chain-id' }],
        size: 16,
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
          pdb: chainsPdb(4, 12),
          chains: 4,
          // sphere impostors (`drawIndirect`) and meshes (`drawIndexedIndirect`)
          // colored by chain: a segment drawn with another one's instances (their segment data) shows
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
