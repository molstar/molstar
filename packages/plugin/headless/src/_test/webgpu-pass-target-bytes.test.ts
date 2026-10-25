/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { describe, expect, it } from 'vitest';
import { now } from '@molstar/core/util/now';
import { MarkerAction } from '@molstar/core/util/marker-action';
import { ReleaseDelayMs } from '@molstar/graphics/canvas3d/passes/pass-targets';
import { HeadlessPluginContext } from '@molstar/plugin-headless/context';
import { DefaultPluginSpec } from '@molstar/plugin/default-spec';
import { WebGPU } from './webgpu-module.js';

const Pdb = `${Array.from({ length: 40 }, (_, i) => {
  const f = (v: number) => v.toFixed(3).padStart(8);
  const serial = String(i + 1).padStart(5);
  const residue = String(i + 1).padStart(4);
  return `ATOM  ${serial}  CA  ALA A${residue}    ${f(Math.cos(i) * 8)}${f(Math.sin(i) * 8)}${f(i * 1.5)}  1.00  0.00           C`;
}).join('\n')}\nEND\n`;

/** every element (`EveryLoci` of the model package, which this package doesn't depend on) */
const EveryLoci = { kind: 'every-loci' as const };

const Size = 64;

const wait = (ms: number) => new Promise((resolve) => setTimeout(resolve, ms));

// the storage of the pass targets per pixel of the drawing buffer, once the targets the props and the scene don't draw
// into are released: a growth here is memory every viewer pays (at 1920×1080, 1 B/px is 2 MiB). The numbers are by
// design; a change that saves memory updates them.
describe.skipIf(!WebGPU.module)('pass target bytes on WebGPU', () => {
  it('stay as planned per scene', async () => {
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
      const representation = await plugin.builders.structure.representation.addRepresentation(structure, {
        type: 'spacefill',
      });
      const canvas = plugin.canvas3d!;
      const passes = plugin.renderer.canvas3dContext!.passes;

      /** pass target bytes per pixel after the frames of the current props and scene and the release */
      const measure = async (iterations = 4) => {
        canvas.commit(true);
        canvas.requestDraw();
        let time = now();
        // the first frame, then the ones accumulating (temporal multi-sampling, tracing) without a change
        for (let i = 0; i < iterations; ++i) {
          canvas.tick((time += 1000), { isSynchronous: true });
          await canvas.graphics.resources.pending.whenReady();
        }
        await wait(ReleaseDelayMs + 200);
        return passes.targets.getByteCount() / (Size * Size);
      };

      const bytes: Record<string, number> = {};
      // defaults: temporal multi-sampling, FXAA, bloom (emissive, nothing emissive), WBOIT
      bytes.opaque = await measure();
      await plugin
        .build()
        .to(representation!)
        .update((old: any) => {
          old.type.params.alpha = 0.5;
        })
        .commit();
      bytes.transparent = await measure();
      plugin.setCanvas3DContextProps({ transparency: 'dpoit' });
      bytes.dpoit = await measure();
      plugin.setCanvas3DContextProps({ transparency: 'wboit' });
      // marking stays once highlighted (hovering is planned)
      canvas.mark({ loci: EveryLoci }, MarkerAction.Highlight);
      bytes.highlighted = await measure();
      canvas.mark({ loci: EveryLoci }, MarkerAction.RemoveHighlight);
      canvas.setProps({ illumination: { ...canvas.props.illumination, enabled: true } });
      bytes.illumination = await measure();
      canvas.setProps({ multiSample: { ...canvas.props.multiSample, mode: 'on' } });
      bytes.illuminationMultiSample = await measure();

      expect(bytes).toEqual({
        // opaque depth-stencil 8, draw color 4, FXAA 4, multi-sample compose (fp16) 8 and first sample 4
        opaque: 28,
        // + WBOIT accumulation (fp16) 8 and weights (r32f) 4, transparency depth-stencil 8
        transparent: 48,
        // WBOIT's in place of: depth peeling's min/max depth (rg32f) 2 × 8, front colors (fp16) 2 × 8, back color 8
        dpoit: 76,
        // + marking layer (fp16) 8, depth, mask and antialiased mask 4 each
        highlighted: 68,
        // the tracing and illumination targets in place of the multi-sampled draw (marking kept)
        illumination: 100,
        // + illumination's multi-sample compose (fp32) 16 and first sample 4
        illuminationMultiSample: 120,
      });
    } finally {
      plugin.dispose();
    }
  }, 60_000);
});
