/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { describe, expect, it } from 'vitest';
import { now } from '@molstar/core/util/now';
import { HiZPass } from '@molstar/graphics/canvas3d/passes/hi-z';
import { createPrepareBuilder } from '@molstar/graphics/canvas3d/passes/pass-builder';
import { Renderer } from '@molstar/graphics/gl/renderer';
import { HeadlessPluginContext } from '@molstar/plugin-headless/context';
import { DefaultPluginSpec } from '@molstar/plugin/default-spec';
import { WebGPU } from './webgpu-module.js';

const Size = 128;

/** a compact grid of atoms in the middle of the view, background around it */
const Pdb = `${Array.from({ length: 48 }, (_, i) => {
  const f = (v: number) => v.toFixed(3).padStart(8);
  const serial = String(i + 1).padStart(5);
  const residue = String(i + 1).padStart(4);
  const [x, y, z] = [(i % 4) * 2.5, (Math.floor(i / 4) % 4) * 2.5, Math.floor(i / 16) * 2.5];
  return `ATOM  ${serial}  CA  ALA A${residue}    ${f(x)}${f(y)}${f(z)}  1.00  0.00           C`;
}).join('\n')}\nEND\n`;

describe.skipIf(!WebGPU.module)('Hi-Z on WebGPU', () => {
  it('levels whose depth is the farthest of the pixels they cover, frame after frame', async () => {
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
      await plugin.builders.structure.representation.addRepresentation(structure, { type: 'spacefill' });
      const canvas = plugin.canvas3d!;
      canvas.setProps({ multiSample: { ...canvas.props.multiSample, mode: 'off' } });
      const { graphics } = canvas;
      const passes = plugin.renderer.canvas3dContext!.passes;
      // a pass of its own on the canvas's draw pass: its levels from 1 (4 pixels a texel at this size)
      const hiZ = new HiZPass(graphics, passes.draw, undefined, { enabled: true, minLevel: 1 });
      hiZ.setViewport(0, 0, Size, Size);

      let time = now();
      // built a few times, prepared (planned without drawing) in between as Canvas3D does: a build after a prepare
      // set the first level's size into the second level's uniforms (`uInvSize` shared), wrong from there on
      const prepare = createPrepareBuilder(graphics, Renderer.create(graphics));
      for (let i = 0; i < 4; ++i) {
        canvas.commit(true);
        canvas.requestDraw();
        canvas.tick((time += 1000) as now.Timestamp, { isSynchronous: true });
        await graphics.resources.pending.whenReady();
        hiZ.plan(prepare, canvas.camera);
        hiZ.render(canvas.camera);
      }
      const image = await graphics.readRenderTarget(passes.draw.drawTarget, 0, 0, Size, Size);
      const pyramid = hiZ.gpuPyramid!;
      expect(pyramid).toBeDefined();
      const { texture } = pyramid;
      const width = texture.getWidth();
      const height = texture.getHeight();
      const levels = new Float32Array(width * height);
      const reader = graphics.createTextureReader({ format: 'r32f' });
      reader.start([texture], 0, 0, width, height);
      for (let t = 0; t < 500 && reader.check([levels]) === 'pending'; ++t) {
        await new Promise((resolve) => setTimeout(resolve, 5));
      }
      reader.dispose();

      // level `minLevel` is the first in the texture: a texel with a depth (not the background's 0) has geometry in all
      // the pixels it covers, the farthest depth of them
      const ts = 1 << (pyramid.levelCount - 1 - pyramid.minLevel);
      const px = Size / ts;
      const background = [image[0], image[1], image[2]];
      const isBackground = (x: number, y: number) => {
        const i = (y * Size + x) * 4;
        return image[i] === background[0] && image[i + 1] === background[1] && image[i + 2] === background[2];
      };
      let filled = 0;
      const wrong: number[][] = [];
      for (let ty = 0; ty < ts; ++ty) {
        for (let tx = 0; tx < ts; ++tx) {
          if (!(levels[ty * width + tx] > 0)) continue;
          filled += 1;
          let covered = true;
          for (let y = 0; y < px && covered; ++y) {
            for (let x = 0; x < px && covered; ++x) if (isBackground(tx * px + x, ty * px + y)) covered = false;
          }
          if (!covered) wrong.push([tx, ty]);
        }
      }
      expect(filled).toBeGreaterThan(4);
      expect(wrong).toEqual([]);
      hiZ.dispose();
    } finally {
      plugin.dispose();
    }
  }, 60_000);
});
