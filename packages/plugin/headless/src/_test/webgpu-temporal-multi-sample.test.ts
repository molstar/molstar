/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { afterEach, describe, expect, it, vi } from 'vitest';
import { now } from '@molstar/core/util/now';
import type { WebGPUContext } from '@molstar/graphics/gl/webgpu/context';
import { HeadlessPluginContext } from '@molstar/plugin-headless/context';
import { DefaultPluginSpec } from '@molstar/plugin/default-spec';
import { WebGPU } from './webgpu-module.js';

// Temporal multi-sampling settles to the image of all its samples in one frame (mode `on`: the same offsets and
// weights), in a frame per sample: the first frame after a change, unjittered, is the first sample (every level's
// first offset is none), not drawn again.

const Pdb = `${Array.from({ length: 16 }, (_, i) => {
  const f = (v: number) => v.toFixed(3).padStart(8);
  const serial = String(i + 1).padStart(5);
  const residue = String(i + 1).padStart(4);
  const [x, y, z] = [(i % 4) * 2.5, Math.floor(i / 4) * 2.5, (i % 3) - 1];
  return `ATOM  ${serial}  CA  ALA A${residue}    ${f(x)}${f(y)}${f(z)}  1.00  0.00           C`;
}).join('\n')}\nEND\n`;

const Size = 48;

afterEach(() => vi.unstubAllGlobals());

async function create(mode: 'on' | 'temporal' | 'off', sampleLevel = 2) {
  const plugin = await HeadlessPluginContext.create(
    { webgpu: WebGPU.module! },
    DefaultPluginSpec(),
    { width: Size, height: Size },
    { backend: 'webgpu' },
  );
  await plugin.init();
  const data = await plugin.builders.data.rawData({ data: Pdb, label: 'test' });
  const trajectory = await plugin.builders.structure.parseTrajectory(data, 'pdb');
  const model = await plugin.builders.structure.createModel(trajectory);
  const structure = await plugin.builders.structure.createStructure(model);
  await plugin.builders.structure.representation.addRepresentation(structure, { type: 'spacefill' });
  const canvas = plugin.canvas3d!;
  canvas.setProps({
    postprocessing: { ...canvas.props.postprocessing, enabled: false },
    multiSample: { ...canvas.props.multiSample, mode, sampleLevel },
  });
  return plugin;
}

/** The frames `plugin`'s canvas draws from a forced draw on, until it draws no more (at most `max`), and its image */
async function settle(plugin: HeadlessPluginContext, max = 20) {
  const canvas = plugin.canvas3d!;
  let time = now();
  const tick = () => canvas.tick((time += 1000) as now.Timestamp, { isSynchronous: true });
  canvas.requestDraw();
  tick();
  await canvas.graphics.resources.pending.whenReady();
  let frames = 0;
  const sub = canvas.didDraw.subscribe(() => (frames += 1));
  canvas.requestDraw();
  for (let i = 0; i < max; ++i) tick();
  sub.unsubscribe();
  // Node has no `ImageData` (restored after the file)
  vi.stubGlobal(
    'ImageData',
    class {
      constructor(
        readonly data: Uint8ClampedArray,
        readonly width: number,
        readonly height: number,
      ) {}
    },
  );
  return { frames, image: await (canvas.graphics as WebGPUContext).readDrawingBuffer() };
}

describe.skipIf(!WebGPU.module)('temporal multi-sampling on WebGPU', () => {
  // level 2: four samples, a frame each; level 5: 32, eight a frame (the running average's weights at their smallest).
  // Mode `on` weights its samples unevenly to spread 8 bit rounding (`roundingRange`): at level 2 by 0.25 ± 0.016, at
  // level 5 from 1/64 to 3/64, the average temporal takes farther from it (measured: hold before the running average
  // as far, max 26, 250 values more than 2 apart, mean 0.36; the running average max 27, 248, mean 0.34)
  for (const [level, samples, perFrame, maxDiff, differShare] of [
    [2, 4, 1, 6, 0.01],
    [5, 32, 8, 32, 0.04],
  ]) {
    it(`settles to the image of all samples at once, level ${level}`, async () => {
      const all = await create('on', level);
      const temporal = await create('temporal', level);
      try {
        const expected = await settle(all);
        const actual = await settle(temporal);
        // the first sample the unjittered frame both draw; the jittered ones `perFrame` a frame (the canvas's own
        // frames, as many in either mode, aside)
        expect(actual.frames - expected.frames).toBe(Math.ceil((samples - 1) / perFrame));
        const a = actual.image.data;
        const e = expected.image.data;
        // the two accumulate apart (the order, the precision, the weights): level 2, 25 of 9216 values more than 2
        // apart, by 4 at most; a sample missing or drawn twice is far more
        let differ = 0;
        let maxDelta = 0;
        for (let i = 0; i < a.length; ++i) {
          const d = Math.abs(a[i] - e[i]);
          if (d > 2) differ += 1;
          maxDelta = Math.max(maxDelta, d);
        }
        expect(maxDelta).toBeLessThanOrEqual(maxDiff);
        expect(differ).toBeLessThanOrEqual(a.length * differShare);
        expect(new Set(e).size).toBeGreaterThan(4);
      } finally {
        all.dispose();
        temporal.dispose();
      }
    });
  }

  // preparing covers what the samples draw with (the running average's constant blend): a draw skipped while its
  // pipeline compiles would lose its sample
  it('settles without skipping a draw once prepared', async () => {
    const temporal = await create('temporal');
    try {
      const canvas = temporal.canvas3d!;
      await canvas.prepare();
      await canvas.graphics.resources.pending.whenReady();
      const skipped = canvas.graphics.stats.skippedDraws;
      let time = now();
      canvas.requestDraw();
      for (let i = 0; i < 6; ++i) canvas.tick((time += 1000) as now.Timestamp, { isSynchronous: true });
      expect(canvas.graphics.stats.skippedDraws).toBe(skipped);
    } finally {
      temporal.dispose();
    }
  });

  // a moving frame (each forced one restarts the samples) presents the draw pass's color target, hold and compose
  // written only when the next sample follows: its image the unjittered scene's, as without multi-sampling
  it('draws a moving frame as without multi-sampling', async () => {
    const off = await create('off');
    const temporal = await create('temporal');
    try {
      const expected = await moving(off);
      const actual = await moving(temporal);
      let maxDelta = 0;
      for (let i = 0; i < actual.data.length; ++i) {
        maxDelta = Math.max(maxDelta, Math.abs(actual.data[i] - expected.data[i]));
      }
      expect(maxDelta).toBe(0);
      expect(new Set(expected.data).size).toBeGreaterThan(4);
    } finally {
      off.dispose();
      temporal.dispose();
    }
  });
});

/** The image of `plugin`'s canvas after forced frames, each a moving one (the last drawn, the first ones compiling) */
async function moving(plugin: HeadlessPluginContext) {
  const canvas = plugin.canvas3d!;
  let time = now();
  for (let i = 0; i < 3; ++i) {
    await canvas.graphics.resources.pending.whenReady();
    canvas.requestDraw();
    canvas.tick((time += 1000) as now.Timestamp, { isSynchronous: true });
  }
  vi.stubGlobal(
    'ImageData',
    class {
      constructor(
        readonly data: Uint8ClampedArray,
        readonly width: number,
        readonly height: number,
      ) {}
    },
  );
  return (canvas.graphics as WebGPUContext).readDrawingBuffer();
}
