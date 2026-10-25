/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { afterAll, beforeAll, describe, expect, it } from 'vitest';
import { now } from '@molstar/core/util/now';
import { ValueCell } from '@molstar/core/util/value-cell';
import * as Shape from '@molstar/graphics/geo/shape/shape';
import { Text } from '@molstar/graphics/geo/geometry/text/text';
import { TextBuilder } from '@molstar/graphics/geo/geometry/text/text-builder';
import { ShapeRepresentation } from '@molstar/graphics/repr/shape/representation';
import { clearFontAtlasCache, setCanvasModule } from '@molstar/graphics/geo/geometry/text/font-atlas';
import type { TextureImage } from '@molstar/graphics/gl/renderable/util';
import { TextSchema } from '@molstar/graphics/gl/renderable/text';
import { createFrame } from '@molstar/graphics/gl/renderable';
import { Renderer } from '@molstar/graphics/gl/renderer';
import { Scene } from '@molstar/graphics/gl/scene';
import { Camera } from '@molstar/graphics/canvas3d/camera';
import { Viewport } from '@molstar/graphics/canvas3d/camera/util';
import type { WebGPUContext } from '@molstar/graphics/gl/webgpu/context';
import { asWebGPUTexture } from '@molstar/graphics/gl/webgpu/texture';
import { HeadlessPluginContext } from '@molstar/plugin-headless/context';
import { DefaultPluginSpec } from '@molstar/plugin/default-spec';
import { WebGPU } from './webgpu-module.js';

// the twenty amino acids (their labels: many glyphs), in chains X and Z (letters none of them has)
const Residues = 'ALA ARG ASN ASP CYS GLN GLU GLY HIS ILE LEU LYS MET PHE PRO SER THR TRP TYR VAL'.split(' ');
const Pdb = `${Residues.map((name, i) => {
  const f = (v: number) => v.toFixed(3).padStart(8);
  const serial = String(i + 1).padStart(5);
  const chain = i < 10 ? 'X' : 'Z';
  const residue = String(i + 1).padStart(4);
  const [x, y, z] = [(i % 10) * 3.8, i < 10 ? 0 : 6, 0];
  return `ATOM  ${serial}  CA  ${name} ${chain}${residue}    ${f(x)}${f(y)}${f(z)}  1.00  0.00           C`;
}).join('\n')}\nEND\n`;

const Size = 64;

// a font whose glyphs differ from each other (no canvas in Node): what the atlas holds is compared, not looked at
beforeAll(() => {
  let glyph = 0;
  const context = {
    // wide glyphs (the atlas caps them): a row of the atlas holds few, the labels' glyphs fill several
    measureText: (text: string) => ({
      width: 1000 * text.length,
      actualBoundingBoxLeft: 0,
      actualBoundingBoxRight: 1000,
    }),
    clearRect() {},
    fillText() {
      glyph += 1;
    },
    getImageData: (_x: number, _y: number, w: number, h: number) => {
      const data = new Uint8ClampedArray(w * h * 4);
      for (let i = 0; i < data.length; ++i) data[i] = (glyph * 37 + i * 11) & 255;
      return { data };
    },
  };
  setCanvasModule({ createCanvas: () => ({ getContext: () => context }) });
});
afterAll(() => {
  setCanvasModule(undefined);
  clearFontAtlasCache();
});

// `GPUBufferUsage` and `GPUMapMode` are not globals in Node
const BufferMapRead = 0x1;
const BufferCopyDst = 0x8;
const MapRead = 0x1;

/** The contents of an `r8unorm` texture, rows of `width` bytes */
async function readR8(graphics: WebGPUContext, texture: GPUTexture, width: number, height: number) {
  const { device } = graphics;
  // rows of a texture-to-buffer copy are aligned to 256 bytes
  const rowBytes = Math.ceil(width / 256) * 256;
  const buffer = device.createBuffer({
    size: rowBytes * height,
    usage: BufferCopyDst | BufferMapRead,
  });
  const encoder = device.createCommandEncoder();
  encoder.copyTextureToBuffer({ texture }, { buffer, bytesPerRow: rowBytes }, { width, height });
  device.queue.submit([encoder.finish()]);
  await buffer.mapAsync(MapRead);
  const padded = new Uint8Array(buffer.getMappedRange());
  const data = new Uint8Array(width * height);
  for (let y = 0; y < height; ++y) data.set(padded.subarray(y * rowBytes, y * rowBytes + width), y * width);
  buffer.unmap();
  buffer.destroy();
  return data;
}

// the texts of one font share its atlas (`.v6/plans/gl/shared-textures.md`): one texture, new glyphs uploaded as the
// rows written, replayed draws of the other texts still bound to it
describe.skipIf(!WebGPU.module)('shared textures on WebGPU', () => {
  it('two label representations share the atlas, new glyphs uploaded as rows', async () => {
    const plugin = await HeadlessPluginContext.create(
      { webgpu: WebGPU.module! },
      DefaultPluginSpec(),
      { width: Size, height: Size },
      { backend: 'webgpu' },
    );
    try {
      await plugin.init();
      const data = await plugin.builders.data.rawData({ data: Pdb, label: 'test' });
      const trajectory = await plugin.builders.structure.parseTrajectory(data, 'pdb');
      const model = await plugin.builders.structure.createModel(trajectory);
      const structure = await plugin.builders.structure.createStructure(model);
      const { representation } = plugin.builders.structure;
      const canvas = plugin.canvas3d!;
      canvas.setProps({
        multiSample: { ...canvas.props.multiSample, mode: 'off' },
        renderer: { ...canvas.props.renderer, drawLists: true },
      });
      const graphics = canvas.graphics as WebGPUContext;
      let time = now();
      const frame = async () => {
        canvas.commit(true);
        canvas.requestDraw();
        canvas.tick((time += 1000) as now.Timestamp, { isSynchronous: true });
        await graphics.resources.pending.whenReady();
      };
      const texts = () => canvas.getRenderObjects().filter((o) => o.type === 'text');

      // the residues' labels: the atlas's first rows of glyphs and more
      await representation.addRepresentation(structure, { type: 'label', typeParams: { level: 'residue' } });
      await frame();
      await frame();
      expect(texts().length).toBe(1);
      const atlas = texts()[0].values.tFont.ref.value as TextureImage<Uint8Array>;
      expect(graphics.sharedTextures.count).toBe(1);
      const held = graphics.sharedTextures.acquire(TextSchema.tFont, atlas);
      const texture = asWebGPUTexture(held.texture).getGPUTexture()!;
      held.release();
      const resources = graphics.drawResourcesVersion;

      // the chains' labels: another text of the same atlas, its new glyphs written in later rows
      const writes: { y: number; rows: number }[] = [];
      const { queue } = graphics.device;
      const writeTexture = queue.writeTexture.bind(queue);
      queue.writeTexture = (destination, data, layout, size) => {
        if (destination.texture === texture) {
          const y = (destination.origin as GPUOrigin3DDict | undefined)?.y ?? 0;
          writes.push({ y, rows: (size as GPUExtent3DDict).height ?? 1 });
        }
        writeTexture(destination, data, layout, size);
      };
      try {
        await representation.addRepresentation(structure, { type: 'label', typeParams: { level: 'chain' } });
        await frame();
        await frame();
      } finally {
        queue.writeTexture = writeTexture;
      }

      expect(texts().length).toBe(2);
      for (const t of texts()) expect(t.values.tFont.ref.value).toBe(atlas);
      expect(graphics.sharedTextures.count).toBe(1);
      // the rows of the new glyphs, not the atlas
      expect(writes.length).toBeGreaterThan(0);
      for (const { rows } of writes) expect(rows).toBeLessThan(atlas.height / 2);
      expect(writes.some(({ y }) => y > 0)).toBe(true);
      // the atlas kept its storage: the first text's recorded draws replay
      expect(asWebGPUTexture(held.texture).getGPUTexture()).toBe(texture);
      expect(graphics.drawResourcesVersion).toBe(resources);
      // what the GPU holds is the atlas
      const uploaded = await readR8(graphics, texture, atlas.width, atlas.height);
      expect(uploaded.length).toBe(atlas.array.length);
      let differing = 0;
      for (let i = 0; i < uploaded.length; ++i) if (uploaded[i] !== atlas.array[i]) differing += 1;
      expect(differing).toBe(0);
      // the atlas the texture's source after the rows (copied from on a backend switch)
      expect(asWebGPUTexture(held.texture).getSource().data).toBe(atlas);
      expect(graphics.takeErrors?.() ?? []).toEqual([]);
    } finally {
      plugin.dispose();
    }
  }, 60_000);

  it('an upload replacing the atlas storage records again what its other holders recorded', async () => {
    const plugin = await HeadlessPluginContext.create(
      { webgpu: WebGPU.module! },
      DefaultPluginSpec(),
      { width: Size, height: Size },
      { backend: 'webgpu' },
    );
    try {
      await plugin.init();
      const canvas = plugin.canvas3d!;
      canvas.setProps({
        multiSample: { ...canvas.props.multiSample, mode: 'off' },
        renderer: { ...canvas.props.renderer, drawLists: true },
      });
      const graphics = canvas.graphics as WebGPUContext;
      let time = now();
      const frame = async () => {
        canvas.commit(true);
        canvas.requestDraw();
        canvas.tick((time += 1000) as now.Timestamp, { isSynchronous: true });
        await graphics.resources.pending.whenReady();
      };

      // more texts than a chunk of a draw list holds (64): a change of one leaves another chunk's recording current
      const props = { sizeFactor: 1 };
      const reprs = [];
      for (let i = 0; i < 65; ++i) {
        const builder = TextBuilder.create(props, 4, 1);
        builder.add(String.fromCharCode(65 + (i % 26)), (i % 9) - 4, Math.floor(i / 9) - 3, 0, 1, 1, 0);
        const text = builder.getText();
        const repr = ShapeRepresentation<{}, Text, Text.Params>(
          () =>
            Shape.create(
              'text',
              {},
              text,
              () => 0xff8800,
              () => 1,
              () => 'text',
            ),
          Text.Utils,
        );
        await repr.createOrUpdate(props, {}).run();
        canvas.add(repr);
        reprs.push(repr);
      }
      canvas.requestCameraReset();
      // recorded, then replayed (drawn directly for the submissions they change in)
      for (let i = 0; i < 6; ++i) await frame();
      const texts = canvas.getRenderObjects().filter((o) => o.type === 'text');
      expect(texts.length).toBe(65);
      const atlas = texts[0].values.tFont.ref.value as TextureImage<Uint8Array>;
      for (const t of texts) expect(t.values.tFont.ref.value).toBe(atlas);
      expect(graphics.sharedTextures.count).toBe(1);
      const { stats } = graphics;
      const stable = { ...stats.drawLists };
      await frame();
      expect(stats.drawLists.replayed).toBeGreaterThan(stable.replayed);
      expect(stats.drawLists.recorded).toBe(stable.recorded);
      expect(graphics.takeErrors?.() ?? []).toEqual([]);

      // the atlas grown in place and one text's cell updated: the upload replaces the storage all texts bind
      const resources = graphics.drawResourcesVersion;
      const grown = new Uint8Array(atlas.array.length * 2);
      grown.set(atlas.array);
      const resized = atlas as { array: Uint8Array; height: number; version?: number };
      resized.array = grown;
      resized.height = atlas.height * 2;
      resized.version = (atlas.version ?? 0) + 1;
      ValueCell.update(texts[0].values.tFont, atlas);
      canvas.update(reprs[0]);
      const before = { ...stats.drawLists };
      await frame();
      await frame();
      expect(graphics.drawResourcesVersion).toBeGreaterThan(resources);
      // nothing replayed what was recorded with the old storage
      expect(stats.drawLists.chunksRecorded + stats.drawLists.chunksDirect).toBeGreaterThan(
        before.chunksRecorded + before.chunksDirect + 1,
      );
      expect(graphics.takeErrors?.() ?? []).toEqual([]);
    } finally {
      plugin.dispose();
    }
  }, 60_000);

  it('a text given another atlas holds its texture, each texture counted once and freed with its last holder', async () => {
    const plugin = await HeadlessPluginContext.create(
      { webgpu: WebGPU.module! },
      DefaultPluginSpec(),
      { width: Size, height: Size },
      { backend: 'webgpu' },
    );
    try {
      await plugin.init();
      const graphics = plugin.canvas3d!.graphics as WebGPUContext;
      /** a text's render object, of the font atlas of `fontQuality` */
      const createText = async (char: string, x: number, fontQuality: number) => {
        const props = { sizeFactor: 1, fontQuality };
        const builder = TextBuilder.create(props, 4, 1);
        builder.add(char, x, 0, 0, 1, 1, 0);
        const text = builder.getText();
        const repr = ShapeRepresentation<{}, Text, Text.Params>(
          () =>
            Shape.create(
              'text',
              {},
              text,
              () => 0xff8800,
              () => 1,
              () => 'text',
            ),
          Text.Utils,
          { graphics },
        );
        await repr.createOrUpdate(props, {}).run();
        return repr.renderObjects[0];
      };
      const a = await createText('A', -2, 2);
      const b = await createText('B', 2, 2);
      // another font quality: another atlas
      const c = await createText('C', 0, 1);
      const atlas = a.values.tFont.ref.value as TextureImage<Uint8Array>;
      const other = c.values.tFont.ref.value as TextureImage<Uint8Array>;
      expect(b.values.tFont.ref.value).toBe(atlas);
      expect(other).not.toBe(atlas);

      const scene = Scene.create(graphics);
      scene.add(a);
      scene.add(b);
      scene.commit();
      const renderer = Renderer.create(graphics);
      const target = graphics.createRenderTarget(Size, Size, 'depth-stencil');
      const camera = new Camera();
      Viewport.set(camera.viewport, 0, 0, Size, Size);
      /** the scene drawn into the target (its pipelines prepared first, none skipped), its draws counted */
      const draw = async () => {
        const { center, radius } = scene.boundingSphere;
        camera.setState(camera.getFocus(center, radius), 0);
        camera.update();
        const descriptor = graphics.describeTarget(target);
        for (const pass of [{ name: 'opaque' }, { name: 'blendedTransparent' }] as const) {
          await graphics.resources.preparePipelines(renderer.collect(pass, scene, descriptor));
        }
        const { stats } = graphics;
        let draws = 0;
        for (let i = 0; i < 2; ++i) {
          const skipped = stats.skippedDraws;
          graphics.beginRenderPass(target);
          renderer.update(camera, scene, createFrame());
          renderer.setViewport(0, 0, Size, Size);
          renderer.renderBlended(scene, camera);
          graphics.flush();
          await graphics.resources.pending.whenReady();
          expect(stats.skippedDraws).toBe(skipped);
          // the frame's calls (reset by `renderer.update`)
          draws = stats.calls.drawInstanced + stats.calls.drawInstancedBase;
        }
        return draws;
      };
      const held = (image: TextureImage<Uint8Array>) => {
        const share = graphics.sharedTextures.acquire(TextSchema.tFont, image);
        share.release();
        return share.texture;
      };
      const textureBytes = (image: TextureImage<Uint8Array>) => image.width * image.height;

      // one texture for both, its bytes split between them
      expect(graphics.sharedTextures.count).toBe(1);
      expect(await draw()).toBeGreaterThanOrEqual(2);
      const [ra, rb] = scene.renderables;
      expect(held(atlas).getByteCount()).toBe(textureBytes(atlas));
      const shareA = ra.getByteCount() + rb.getByteCount();
      expect(shareA).toBeGreaterThanOrEqual(textureBytes(atlas));
      expect(shareA).toBeLessThan(2 * textureBytes(atlas));
      expect(graphics.takeErrors?.() ?? []).toEqual([]);

      // the second text given the other atlas: it holds the texture of that one, the first keeps its own
      ValueCell.update(b.values.tFont, other);
      scene.update([b], false);
      expect(graphics.sharedTextures.count).toBe(2);
      expect(await draw()).toBeGreaterThanOrEqual(2);
      expect(rb.getByteCount()).toBeGreaterThanOrEqual(textureBytes(other));
      expect(rb.getByteCount()).toBeLessThan(textureBytes(atlas));
      expect(ra.getByteCount()).toBeGreaterThanOrEqual(textureBytes(atlas));
      expect(held(other)).not.toBe(held(atlas));
      expect(graphics.takeErrors?.() ?? []).toEqual([]);

      // removed, each texture freed with its last holder; the remaining text still drawn
      const textures = graphics.stats.resourceCounts.texture;
      scene.remove(a);
      scene.commit();
      expect(graphics.sharedTextures.count).toBe(1);
      expect(graphics.stats.resourceCounts.texture).toBeLessThan(textures);
      expect(await draw()).toBeGreaterThanOrEqual(1);
      scene.remove(b);
      scene.commit();
      expect(graphics.sharedTextures.count).toBe(0);
      expect(graphics.takeErrors?.() ?? []).toEqual([]);

      renderer.dispose();
      target.destroy();
      scene.clear();
    } finally {
      plugin.dispose();
    }
  }, 60_000);
});
