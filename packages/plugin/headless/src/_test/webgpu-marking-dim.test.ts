/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { describe, expect, it, vi } from 'vitest';
import { now } from '@molstar/core/util/now';
import { MarkerAction } from '@molstar/core/util/marker-action';
import { HeadlessPluginContext } from '@molstar/plugin-headless/context';
import { DefaultPluginSpec } from '@molstar/plugin/default-spec';
import type { WebGPUContext } from '@molstar/graphics/gl/webgpu/context';
import { WebGPU } from './webgpu-module.js';

/** two chains of a grid of atoms, side by side: an object each, chain A on the left */
const TwoChainsPdb = `${Array.from({ length: 48 }, (_, i) => {
  const f = (v: number) => v.toFixed(3).padStart(8);
  const serial = String(i + 1).padStart(5);
  const chain = i < 24 ? 'A' : 'B';
  const residue = String((i % 24) + 1).padStart(4);
  const [x, y, z] = [(i % 4) * 2.5 + (i < 24 ? 0 : 12), (Math.floor(i / 4) % 6) * 2.5, 0];
  return `ATOM  ${serial}  CA  ALA ${chain}${residue}    ${f(x)}${f(y)}${f(z)}  1.00  0.00           C`;
}).join('\n')}\nEND\n`;

const Size = 64;

describe.skipIf(!WebGPU.module)('dimming by the marking pass on WebGPU', () => {
  // dimming set after the canvas was created reached the marking pass (the frames' props kept the renderer's as created),
  // and transparent unmarked objects are dimmed (the marking depth taken from the scene's opaque depth lacked them)
  it('dims the unmarked transparent chain, the dimming set on the created canvas', async () => {
    const plugin = await HeadlessPluginContext.create(
      { webgpu: WebGPU.module! },
      DefaultPluginSpec(),
      { width: Size, height: Size },
      { backend: 'webgpu' },
    );
    // Node has no `ImageData` (`readDrawingBuffer` returns one); the setup restores stubbed globals after the file
    if (typeof ImageData === 'undefined') {
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
    }
    try {
      await plugin.init();
      const data = await plugin.builders.data.rawData({ data: TwoChainsPdb, label: 'test' });
      const trajectory = await plugin.builders.structure.parseTrajectory(data, 'pdb');
      const model = await plugin.builders.structure.createModel(trajectory);
      const structure = await plugin.builders.structure.createStructure(model);
      // as the default preset's water: transparent, less than the marking's threshold for hiding
      await plugin.builders.structure.representation.addRepresentation(structure, {
        type: 'spacefill',
        typeParams: { alpha: 0.6 },
      });
      const canvas = plugin.canvas3d!;
      canvas.setProps({
        multiSample: { ...canvas.props.multiSample, mode: 'off' },
        marking: { ...canvas.props.marking, enabled: true },
        renderer: { ...canvas.props.renderer, colorMarker: true, dimStrength: 0 },
      });
      let time = now();
      const frame = async () => {
        canvas.requestDraw();
        canvas.tick((time += 1000) as now.Timestamp, { isSynchronous: true });
        await canvas.graphics.resources.pending.whenReady();
        return (await (canvas.graphics as WebGPUContext).readDrawingBuffer()).data;
      };
      canvas.commit(true);
      await frame();
      const objects = canvas.getRenderObjects().filter((o) => o.state.visible);
      expect(objects.length).toBe(2);
      // an atom of chain A selected: chain B unmarked
      canvas.mark(canvas.getLoci({ objectId: objects[0].id, instanceId: 0, groupId: 0 }), MarkerAction.Select);
      for (let i = 0; i < 3; ++i) await frame();
      const undimmed = await frame();

      canvas.setProps({ renderer: { ...canvas.props.renderer, dimColor: 0xffffff, dimStrength: 0.8 } });
      for (let i = 0; i < 3; ++i) await frame();
      const dimmed = await frame();

      /** the mean luminance of the right half (chain B) */
      const right = (pixels: Uint8ClampedArray) => {
        let sum = 0;
        let count = 0;
        for (let y = 0; y < Size; ++y) {
          for (let x = Size / 2; x < Size; ++x) {
            const o = (y * Size + x) * 4;
            sum += 0.299 * pixels[o] + 0.587 * pixels[o + 1] + 0.114 * pixels[o + 2];
            count += 1;
          }
        }
        return sum / count;
      };
      // toward the white dim color
      expect(right(dimmed)).toBeGreaterThan(right(undimmed) + 10);
    } finally {
      plugin.dispose();
    }
  }, 60_000);
});
