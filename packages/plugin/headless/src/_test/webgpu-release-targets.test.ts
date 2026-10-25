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

const wait = (ms: number) => new Promise((resolve) => setTimeout(resolve, ms));

// pass targets take memory while the props and the scene draw into them: a pass turned off releases its targets a
// while after the change, those of the passes still drawn (and of hovering: marking) stay; an image pass releases its
// targets a while after its last render
describe.skipIf(!WebGPU.module)('releasing pass targets on WebGPU', () => {
  it('releases the targets of passes turned off, and of image passes after rendering', async () => {
    const plugin = await HeadlessPluginContext.create(
      { webgpu: WebGPU.module! },
      DefaultPluginSpec(),
      { width: 64, height: 64 },
      { backend: 'webgpu' },
    );
    try {
      await plugin.init();
      const data = await plugin.builders.data.rawData({ data: Pdb, label: 'test' });
      const trajectory = await plugin.builders.structure.parseTrajectory(data, 'pdb');
      const model = await plugin.builders.structure.createModel(trajectory);
      const structure = await plugin.builders.structure.createStructure(model);
      await plugin.builders.structure.representation.addRepresentation(structure, {
        type: 'spacefill',
        typeParams: { alpha: 0.5 },
      });
      const canvas = plugin.canvas3d!;
      const passes = plugin.renderer.canvas3dContext!.passes;
      const frames = async () => {
        for (let i = 0; i < 2; ++i) {
          canvas.commit(true);
          canvas.requestDraw();
          canvas.tick(now(), { isSynchronous: true });
          await canvas.graphics.resources.pending.whenReady();
        }
      };

      await frames();
      // transparent: drawn with WBOIT
      expect(passes.draw.wboit.getByteCount()).toBeGreaterThan(0);
      expect(passes.illumination.getByteCount()).toBe(0);

      canvas.setProps({ illumination: { ...canvas.props.illumination, enabled: true } });
      await frames();
      expect(passes.illumination.getByteCount()).toBeGreaterThan(0);

      canvas.setProps({ illumination: { ...canvas.props.illumination, enabled: false } });
      await frames();
      // not yet: released a while after the change
      expect(passes.illumination.getByteCount()).toBeGreaterThan(0);
      await wait(ReleaseDelayMs + 200);
      expect(passes.illumination.getByteCount()).toBe(0);
      // still drawn into
      expect(passes.draw.wboit.getByteCount()).toBeGreaterThan(0);
      // marked, then not: the marking targets stay, kept for the next hover
      canvas.mark({ loci: EveryLoci }, MarkerAction.Highlight);
      await frames();
      const marking = passes.draw.marking.getByteCount();
      expect(marking).toBeGreaterThan(0);
      canvas.mark({ loci: EveryLoci }, MarkerAction.RemoveHighlight);
      canvas.setProps({ illumination: { ...canvas.props.illumination, enabled: false } });
      await frames();
      await wait(ReleaseDelayMs + 200);
      expect(passes.draw.marking.getByteCount()).toBe(marking);

      // an image pass releases its targets a while after its last render
      await plugin.getImageRaw({ width: 64, height: 64 });
      const { imagePass } = plugin.renderer;
      expect(imagePass.getByteCount()).toBeGreaterThan(0);
      await wait(ReleaseDelayMs + 200);
      expect(imagePass.getByteCount()).toBe(0);
      // and allocates them again
      await plugin.getImageRaw({ width: 64, height: 64 });
      expect(imagePass.getByteCount()).toBeGreaterThan(0);
    } finally {
      plugin.dispose();
    }
  });
});
