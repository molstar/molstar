/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { describe, expect, it, vi } from 'vitest';
import { now } from '@molstar/core/util/now';
import { Vec3 } from '@molstar/core/math/linear-algebra/3d/vec3';
import { HiZPass } from '@molstar/graphics/canvas3d/passes/hi-z';
import { HeadlessPluginContext } from '@molstar/plugin-headless/context';
import { DefaultPluginSpec } from '@molstar/plugin/default-spec';
import type { RenderTarget } from '@molstar/graphics/gl/render-target';
import type { WebGPUContext } from '@molstar/graphics/gl/webgpu/context';
import { WebGPU } from './webgpu-module.js';

const Size = 64;

/**
 * A solid wall of atoms (chain A, 30x30, at z = 0) and a small cube of atoms (chain B) as 20 copies in a row behind its
 * middle along -z: a merged object culled on the GPU, the copies hidden by the wall from the front, in view from the side
 */
function wallPdb() {
  const f = (v: number) => v.toFixed(3).padStart(8);
  const atoms: string[] = [];
  const atom = (chain: string, j: number, x: number, y: number, z: number) =>
    atoms.push(
      `ATOM  ${String(atoms.length + 1).padStart(5)}  CA  ALA ${chain}${String(j + 1).padStart(4)}    ${f(x)}${f(y)}${f(z)}  1.00  0.00           C`,
    );
  for (let j = 0; j < 900; ++j) atom('A', j, ((j % 30) - 14.5) * 1.2, (Math.floor(j / 30) - 14.5) * 1.2, 0);
  for (let j = 0; j < 8; ++j)
    atom('B', j, ((j % 2) - 0.5) * 1.2, ((Math.floor(j / 2) % 2) - 0.5) * 1.2, -Math.floor(j / 4) * 1.2);
  const identity = (n: string, z: number) => [
    `REMARK 350   BIOMT1 ${n}  1.000000  0.000000  0.000000        0.00000`,
    `REMARK 350   BIOMT2 ${n}  0.000000  1.000000  0.000000        0.00000`,
    `REMARK 350   BIOMT3 ${n}  0.000000  0.000000  1.000000 ${f(z)}`,
  ];
  const biomt: string[] = ['REMARK 350 BIOMOLECULE: 1', 'REMARK 350 APPLY THE FOLLOWING TO CHAINS: A'];
  biomt.push(...identity('   1', 0));
  biomt.push('REMARK 350 APPLY THE FOLLOWING TO CHAINS: B');
  for (let k = 0; k < 20; ++k) biomt.push(...identity(String(k + 2).padStart(4), -10 - 4 * k));
  return `${biomt.join('\n')}\n${atoms.join('\n')}\nEND\n`;
}

describe.skipIf(!WebGPU.module)('illumination inputs culled by occlusion on WebGPU', () => {
  // the inputs (tracing, back depth) culled on the GPU against the previous frame's pyramid, and in a second phase
  // against that of the tracing input's depth so far: a frame after a jump of the camera (the copies behind the front
  // one revealed, the previous pyramid hiding them) the same as without occlusion
  it('draws what is drawn without occlusion, the camera jumping', async () => {
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
    const planFrame = vi.spyOn(HiZPass.prototype, 'planFrame');
    try {
      await plugin.init();
      const data = await plugin.builders.data.rawData({ data: wallPdb(), label: 'test' });
      const trajectory = await plugin.builders.structure.parseTrajectory(data, 'pdb');
      const model = await plugin.builders.structure.createModel(trajectory);
      const structure = await plugin.builders.structure.createStructure(model, {
        name: 'assembly',
        params: { id: '1' },
      });
      await plugin.builders.structure.representation.addRepresentation(structure, {
        type: 'spacefill',
        // a cell per copy and chain
        typeParams: { renderMerged: true, gpuCulling: true, cellSize: 3, batchSize: 12 },
      });
      const canvas = plugin.canvas3d!;
      canvas.setProps({
        multiSample: { ...canvas.props.multiSample, mode: 'off' },
        camera: { ...canvas.props.camera, mode: 'perspective' },
        hiZ: { ...canvas.props.hiZ, enabled: false, gpu: true, twoPhase: true },
        // the quality not adapted to the frame times: the first iteration's fixed, the following ones climbing by a
        // step each (no target frame rate to fall below)
        illumination: { ...canvas.props.illumination, enabled: true, rendersPerFrame: [1, 1], targetFps: 0 },
      });
      // the viewports synced (the Hi-Z pass's among them)
      canvas.handleResize();
      const graphics = canvas.graphics as WebGPUContext;
      const { illumination } = plugin.renderer.canvas3dContext!.passes;
      let time = now();
      /** `count` ticks, forced or not (not: the illumination iterating) */
      const frames = async (count: number, force = true) => {
        for (let i = 0; i < count; ++i) {
          canvas.commit(true);
          if (force) canvas.requestDraw();
          canvas.tick((time += 1000) as now.Timestamp, { isSynchronous: true });
          await graphics.resources.pending.whenReady();
        }
        return Array.from((await graphics.readDrawingBuffer()).data);
      };
      await frames(1);
      const { center, radius } = canvas.boundingSphere;
      const target = Vec3.clone(center);
      const up = Vec3.create(0, 1, 0);
      // close to the wall, which fills the view (the copies behind it); then from the side, in one step
      const front = Vec3.create(center[0], center[1], 15);
      const side = Vec3.add(Vec3(), center, Vec3.create(radius, 0, 0));
      // the back depth's target (private: the trace reads it, nothing shows it)
      const thickness = (illumination as unknown as { tracing: { thicknessTarget: RenderTarget } }).tracing
        .thicknessTarget;
      /**
       * the back depth of the front view (the copies hidden by the wall, their back faces the farthest); the frame after
       * the jump, its back depth; the image of the iterations after it
       */
      const jump = async () => {
        canvas.camera.setState({ target, position: front, up }, 0);
        await frames(4);
        const hidden = Array.from(await graphics.readRenderTarget(thickness, 0, 0, Size, Size));
        illumination.restart(true);
        canvas.camera.setState({ target, position: side, up }, 0);
        const moving = await frames(1);
        const back = Array.from(await graphics.readRenderTarget(thickness, 0, 0, Size, Size));
        return { hidden, moving, back, still: await frames(6, false) };
      };

      planFrame.mockClear();
      const occluded = await jump();
      const secondPhase = planFrame.mock.results.some((r) => r.value !== undefined);

      canvas.setProps({ hiZ: { ...canvas.props.hiZ, gpu: false } });
      const all = await jump();
      expect(new Set(all.moving).size).toBeGreaterThan(2);
      expect(occluded.moving).toEqual(all.moving);
      expect(occluded.back).toEqual(all.back);
      // the back depth not culled by the front's occlusion: what the wall hides is in it
      expect(occluded.hidden).toEqual(all.hidden);
      expect(occluded.still).toEqual(all.still);
      // the second phase ran, against a pyramid of the tracing input's depth
      expect(secondPhase).toBe(true);
    } finally {
      plugin.dispose();
    }
  }, 60_000);
});
