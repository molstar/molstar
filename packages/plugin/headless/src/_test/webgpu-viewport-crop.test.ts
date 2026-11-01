/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { describe, expect, it, vi } from 'vitest';
import { Vec2, Vec3, Vec4 } from '@molstar/core/math/linear-algebra';
import { now } from '@molstar/core/util/now';
import type { WebGPUContext } from '@molstar/graphics/gl/webgpu/context';
import type { Camera } from '@molstar/graphics/canvas3d/camera';
import { HeadlessPluginContext } from '@molstar/plugin-headless/context';
import { DefaultPluginSpec } from '@molstar/plugin/default-spec';
import { WebGPU } from './webgpu-module.js';

const Pdb = `${Array.from({ length: 16 }, (_, i) => {
  const f = (v: number) => v.toFixed(3).padStart(8);
  const serial = String(i + 1).padStart(5);
  const residue = String(i + 1).padStart(4);
  const [x, y, z] = [(i % 4) * 2.5, Math.floor(i / 4) * 2.5, (i % 3) - 1];
  return `ATOM  ${serial}  CA  ALA A${residue}    ${f(x)}${f(y)}${f(z)}  1.00  0.00           C`;
}).join('\n')}\nEND\n`;

const Size = 48;
const Frame = { width: 64, height: 40 };

/** A spacefill scene on a canvas of `width` × `height`, drawn through `viewport` (the canvas where not given) */
async function create(width: number, height: number) {
  const plugin = await HeadlessPluginContext.create(
    { webgpu: WebGPU.module! },
    DefaultPluginSpec(),
    { width, height },
    { backend: 'webgpu' },
  );
  await plugin.init();
  const data = await plugin.builders.data.rawData({ data: Pdb, label: 'test' });
  const trajectory = await plugin.builders.structure.parseTrajectory(data, 'pdb');
  const model = await plugin.builders.structure.createModel(trajectory);
  const structure = await plugin.builders.structure.createStructure(model);
  await plugin.builders.structure.representation.addRepresentation(structure, { type: 'spacefill' });
  const canvas = plugin.canvas3d!;
  // what a pixel shows depends on its own fragments only: no screen-space effects or sample accumulation
  canvas.setProps({
    postprocessing: { ...canvas.props.postprocessing, enabled: false },
    multiSample: { ...canvas.props.multiSample, mode: 'off' },
  });
  return plugin;
}

/** Draws `plugin`'s canvas at `snapshot` once its pipelines are ready; its drawing buffer, rows top to bottom */
async function draw(plugin: HeadlessPluginContext, snapshot: Camera.Snapshot) {
  const canvas = plugin.canvas3d!;
  canvas.camera.setState(snapshot, 0);
  let time = now();
  const frame = () => {
    canvas.requestDraw();
    canvas.tick((time += 1000) as now.Timestamp, { isSynchronous: true });
  };
  frame();
  await canvas.graphics.resources.pending.whenReady();
  frame();
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
  return (canvas.graphics as WebGPUContext).readDrawingBuffer();
}

// a Canvas3D frame extending past the canvas is cropped by Canvas3D (the camera's crop), so that no viewport outside
// the target reaches the backend (Safari's WebGPU rejects one): it draws what the whole frame draws there
describe.skipIf(!WebGPU.module)('a Canvas3D frame extending past the canvas on WebGPU', () => {
  const cases: { name: string; x: number; y: number; width?: number; height?: number }[] = [
    { name: 'off the left and top', x: -16, y: -8 },
    { name: 'off the right and bottom', x: 8, y: 16 },
    { name: 'larger than the canvas, off every side', x: -8, y: -4, width: 64, height: 56 },
  ];
  for (const mode of ['perspective', 'orthographic'] as const)
    for (const c of cases) {
      it(`draws the part inside as the whole frame draws it (${c.name}, ${mode})`, async () => {
        const frame = { width: c.width ?? Frame.width, height: c.height ?? Frame.height };
        const whole = await create(frame.width, frame.height);
        const cropped = await create(Size, Size);
        try {
          const snapshot = { ...whole.canvas3d!.camera.getSnapshot(), mode };
          const expected = await draw(whole, snapshot);
          const canvas = cropped.canvas3d!;
          const pixelRatio = canvas.graphics.pixelRatio;
          canvas.setProps({
            viewport: {
              name: 'static-frame',
              params: {
                x: c.x / pixelRatio,
                y: c.y / pixelRatio,
                width: frame.width / pixelRatio,
                height: frame.height / pixelRatio,
              },
            },
          });
          const actual = await draw(cropped, snapshot);
          expect(canvas.camera.crop.enabled).toBe(true);

          // the frame's pixels inside the canvas
          const x0 = Math.max(c.x, 0);
          const y0 = Math.max(c.y, 0);
          const x1 = Math.min(c.x + frame.width, Size);
          const y1 = Math.min(c.y + frame.height, Size);
          let drawn = 0;
          let differing = 0;
          for (let y = y0; y < y1; ++y) {
            for (let x = x0; x < x1; ++x) {
              const i = (y * Size + x) * 4;
              const j = ((y - c.y) * frame.width + (x - c.x)) * 4;
              let d = 0;
              for (let k = 0; k < 4; ++k) d = Math.max(d, Math.abs(actual.data[i + k] - expected.data[j + k]));
              if (d > 2) differing += 1;
              if (expected.data[j] !== expected.data[(0 * frame.width + 0) * 4]) drawn += 1;
            }
          }
          // the scene reaches into the cropped part, and is drawn there as in the whole frame
          expect(drawn).toBeGreaterThan(50);
          expect(differing).toBe(0);

          // picking at an atom drawn in the cropped part finds it
          const p = canvas.camera.project(Vec4(), Vec3.create(2.5, 2.5, 1));
          expect(p[0]).toBeGreaterThan(x0);
          expect(p[0]).toBeLessThan(x1);
          expect(p[1]).toBeGreaterThan(y0);
          expect(p[1]).toBeLessThan(y1);
          const pick = await canvas.asyncIdentify(Vec2.create(p[0] / pixelRatio, p[1] / pixelRatio))?.wait();
          expect(pick?.id.objectId).toBeGreaterThanOrEqual(0);
        } finally {
          whole.dispose();
          cropped.dispose();
        }
      }, 60_000);
    }
});
