/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { describe, expect, it } from 'vitest';
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

/** chains of 12 atoms side by side: a unit each, merged */
const ChainsPdb = `${Array.from({ length: 48 }, (_, i) => {
  const f = (v: number) => v.toFixed(3).padStart(8);
  const serial = String(i + 1).padStart(5);
  const chain = 'ABCD'[Math.floor(i / 12)];
  const residue = String((i % 12) + 1).padStart(4);
  const [x, y, z] = [(i % 3) * 3.8 + Math.floor(i / 12) * 14, Math.floor((i % 12) / 3) * 3.8, (i % 2) * 1.5];
  return `ATOM  ${serial}  CA  ALA ${chain}${residue}    ${f(x)}${f(y)}${f(z)}  1.00  0.00           C`;
}).join('\n')}\nEND\n`;

const Size = 64;

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
    list.counts[i] = o.members[i].values.drawCount.ref.value;
    list.offsets[i] = segments.elementOffsets[i];
    list.instanceCounts[i] = segments.instanceCounts[i];
    list.baseInstances[i] = segments.instanceBases[i];
  }
  return list;
}

describe.skipIf(!WebGPU.module)('indirect draw lists on WebGPU', () => {
  it('draw what the same direct lists draw, indexed or not, from instances other than the first', async () => {
    const plugin = await HeadlessPluginContext.create(
      { webgpu: WebGPU.module! },
      DefaultPluginSpec(),
      { width: Size, height: Size },
      { backend: 'webgpu' },
    );
    try {
      await plugin.init();
      const graphics = plugin.canvas3d!.graphics;
      expect(graphics.capabilities.indirectInstanceRanges).toBe(true);
      const data = await plugin.builders.data.rawData({ data: ChainsPdb, label: 'test' });
      const trajectory = await plugin.builders.structure.parseTrajectory(data, 'pdb');
      const model = await plugin.builders.structure.createModel(trajectory);
      const structure = await plugin.builders.structure.createStructure(model);
      const { representation } = plugin.builders.structure;
      // sphere impostors (`drawIndirect`) and meshes (`drawIndexedIndirect`)
      // colored by chain: a segment drawn with another one's instances (their segment data) shows
      await representation.addRepresentation(structure, {
        type: 'spacefill',
        typeParams: { renderMerged: true },
        color: 'chain-id',
      });
      await representation.addRepresentation(structure, {
        type: 'ball-and-stick',
        typeParams: { renderMerged: true, sizeFactor: 0.2, tryUseImpostor: false },
        color: 'chain-id',
      });
      plugin.canvas3d!.commit(true);

      const objects = plugin.canvas3d!.getRenderObjects().filter((o) => isMergedRenderObject(o));
      expect(objects.map((o) => o.type).sort()).toEqual(['mesh', 'spheres']);
      const scene = Scene.create(graphics);
      for (const o of objects) scene.add(o);
      scene.commit();
      const camera = new Camera();
      Viewport.set(camera.viewport, 0, 0, Size, Size);
      const { center, radius } = scene.boundingSphere;
      camera.setState(camera.getFocus(center, radius), 0);
      camera.update();

      const renderer = Renderer.create(graphics);
      const target = graphics.createRenderTarget(Size, Size, 'depth-stencil');
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
        expect(list.count).toBe(4);
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
        renderer.setViewport(0, 0, Size, Size);
        renderer.renderOpaque(scene, camera);
        graphics.flush();
        return Array.from(await graphics.readRenderTarget(target, 0, 0, Size, Size));
      };
      // the first frame differs (what it draws with settles in it): from the second on, frames are the same
      await frame();
      const own = await frame();
      const directImage = await frame(direct);
      const indirectImage = await frame(indirect);
      // the objects are drawn (not only the background), the lists draw what the objects' own draw
      expect(new Set(own).size).toBeGreaterThan(2);
      expect(directImage).toEqual(own);
      expect(indirectImage).toEqual(own);
    } finally {
      plugin.dispose();
    }
  }, 60_000);
});
