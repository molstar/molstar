/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { describe, expect, it } from 'vitest';
import { now } from '@molstar/core/util/now';
import { testBundleFragDepth } from '@molstar/graphics/gl/webgpu/bundle-frag-depth';
import { HeadlessPluginContext } from '@molstar/plugin-headless/context';
import { DefaultPluginSpec } from '@molstar/plugin/default-spec';
import { WebGPU } from './webgpu-module.js';

const Pdb = `${Array.from({ length: 16 }, (_, i) => {
  const f = (v: number) => v.toFixed(3).padStart(8);
  const serial = String(i + 1).padStart(5);
  const residue = String(i + 1).padStart(4);
  const [x, y, z] = [(i % 4) * 2.5, Math.floor(i / 4) * 2.5, 0];
  return `ATOM  ${serial}  CA  ALA A${residue}    ${f(x)}${f(y)}${f(z)}  1.00  0.00           C`;
}).join('\n')}\nEND\n`;

const Size = 64;

/** The chunks replayed by 4 frames of spacefill (sphere impostors, which write `frag_depth`) and the last image */
async function render(fragDepthInBundles?: boolean) {
  const plugin = await HeadlessPluginContext.create(
    { webgpu: WebGPU.module! },
    DefaultPluginSpec(),
    { width: Size, height: Size },
    { backend: 'webgpu', webgpu: { fragDepthInBundles } },
  );
  try {
    await plugin.init();
    const data = await plugin.builders.data.rawData({ data: Pdb, label: 'test' });
    const trajectory = await plugin.builders.structure.parseTrajectory(data, 'pdb');
    const model = await plugin.builders.structure.createModel(trajectory);
    const structure = await plugin.builders.structure.createStructure(model);
    await plugin.builders.structure.representation.addRepresentation(structure, { type: 'spacefill' });
    const canvas = plugin.canvas3d!;
    canvas.setProps({ renderer: { ...canvas.props.renderer, drawLists: true } });
    const { stats } = canvas.graphics;
    let time = now();
    const frame = () => {
      canvas.requestDraw();
      canvas.tick((time += 1000) as now.Timestamp, { isSynchronous: true });
    };
    // pipelines (and the probe) settle, then frames record and replay
    frame();
    await canvas.graphics.resources.pending.whenReady();
    const replayed = stats.drawLists.replayed;
    for (let i = 0; i < 4; ++i) frame();
    const image = await plugin.getImageRaw({ width: Size, height: Size });
    return { replayed: stats.drawLists.replayed - replayed, image: Array.from(image.data) };
  } finally {
    plugin.dispose();
  }
}

// Safari loses the depth a fragment writes in draws replayed from a render bundle (after another pipeline, two or more
// uniform buffers at dynamic offsets): such draws are drawn directly there
describe.skipIf(!WebGPU.module)('frag_depth in render bundles on WebGPU', () => {
  it('the probe finds the depth kept on Dawn', async () => {
    // `gpu` kept alive while its device is in use (collecting it crashes Dawn's Node.js bindings)
    const gpu = WebGPU.module!.create([]);
    const adapter = await gpu.requestAdapter();
    const device = await adapter!.requestDevice();
    try {
      expect(await testBundleFragDepth(device)).toBe(true);
    } finally {
      device.destroy();
      expect(gpu).toBeDefined();
    }
  });

  it('draws writing frag_depth replay where the probe keeps the depth, are drawn directly where not', async () => {
    const probed = await render();
    const kept = await render(true);
    const lost = await render(false);
    expect(probed.replayed).toBe(kept.replayed);
    // the spheres' chunks drawn directly in the frames where the bundles lose their depth, replayed where they keep it
    // (how many: the passes drawing them, not counted here)
    expect(lost.replayed).toBeLessThan(kept.replayed);
    expect(lost.image).toEqual(kept.image);
    expect(probed.image).toEqual(kept.image);
  }, 60_000);
});
