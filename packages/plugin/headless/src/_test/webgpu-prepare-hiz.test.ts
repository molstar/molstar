/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { readFileSync } from 'node:fs';
import { describe, expect, it } from 'vitest';
import { now } from '@molstar/core/util/now';
import { HeadlessPluginContext } from '@molstar/plugin-headless/context';
import { DefaultPluginSpec } from '@molstar/plugin/default-spec';
import { WebGPU } from './webgpu-module.js';

const readExample = (name: string) =>
  readFileSync(new URL(`../../../../../data/examples/${name}`, import.meta.url), 'utf8');

// Preparing props a canvas is not set to (`Canvas3D.prepare({ props })`) plans their Hi-Z: on a headless canvas, which
// has its viewport from the start (no resize), the Hi-Z pyramid's pipeline is compiled by preparing, not by the first
// frame drawn with those props
describe.skipIf(!WebGPU.module)('Hi-Z prepared with given props on WebGPU', () => {
  it('compile the Hi-Z pyramid when preparing, nothing of it when drawing', async () => {
    const plugin = await HeadlessPluginContext.create(
      { webgpu: WebGPU.module! },
      DefaultPluginSpec(),
      { width: 64, height: 64 },
      { backend: 'webgpu' },
    );
    try {
      await plugin.init();
      const data = await plugin.builders.data.rawData({ data: readExample('1cbs_updated.cif'), label: '1cbs' });
      const trajectory = await plugin.builders.structure.parseTrajectory(data, 'mmcif');
      const model = await plugin.builders.structure.createModel(trajectory);
      const structure = await plugin.builders.structure.createStructure(model);
      await plugin.builders.structure.representation.addRepresentation(structure, { type: 'spacefill' });

      const canvas = plugin.canvas3d!;
      const { pending } = canvas.graphics.resources;
      expect(canvas.props.hiZ.enabled).toBe(false);
      const hiZ = { ...canvas.props.hiZ, enabled: true };
      const isHiZPipeline = (label: string) => label.startsWith('pipeline hi-z');

      const beforePrepare = pending.trackedCount;
      await canvas.prepare({ props: { hiZ } });
      await pending.whenReady();
      expect(pending.trackedSince(beforePrepare).some(isHiZPipeline)).toBe(true);

      const beforeDraw = pending.trackedCount;
      canvas.setProps({ hiZ });
      for (let i = 0; i < 3; ++i) {
        canvas.requestDraw();
        canvas.tick(now(), { isSynchronous: true });
      }
      await pending.whenReady();
      expect(pending.trackedSince(beforeDraw).filter(isHiZPipeline)).toEqual([]);
    } finally {
      plugin.dispose();
    }
  }, 60_000);
});
