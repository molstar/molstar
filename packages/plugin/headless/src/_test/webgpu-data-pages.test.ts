/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { describe, expect, it, vi } from 'vitest';
import { HeadlessPluginContext } from '@molstar/plugin-headless/context';
import { DefaultPluginSpec } from '@molstar/plugin/default-spec';
import { WebGPU } from './webgpu-module.js';

const Pdb = `ATOM      1  N   ALA A   1      -1.458   0.000   0.000  1.00  0.00           N
ATOM      2  CA  ALA A   1       0.000   0.000   0.000  1.00  0.00           C
ATOM      3  C   ALA A   1       0.551   1.420   0.000  1.00  0.00           C
ATOM      4  O   ALA A   1      -0.200   2.390   0.000  1.00  0.00           O
ATOM      5  CB  ALA A   1       0.540  -0.760  -1.210  1.00  0.00           C
ATOM      6  N   GLY A   2       1.870   1.550   0.000  1.00  0.00           N
ATOM      7  CA  GLY A   2       2.520   2.860   0.000  1.00  0.00           C
ATOM      8  C   GLY A   2       4.030   2.720   0.000  1.00  0.00           C
ATOM      9  O   GLY A   2       4.570   1.610   0.000  1.00  0.00           O
END
`;

const Size = 64;

/** The image of the structure as spacefill, colored by element index, the data binding limited to `maxDataBindingSize` */
async function render(maxDataBindingSize?: number) {
  const plugin = await HeadlessPluginContext.create(
    { webgpu: WebGPU.module! },
    DefaultPluginSpec(),
    { width: Size, height: Size },
    { backend: 'webgpu', webgpu: { maxDataBindingSize } },
  );
  try {
    await plugin.init();
    const data = await plugin.builders.data.rawData({ data: Pdb, label: 'test' });
    const trajectory = await plugin.builders.structure.parseTrajectory(data, 'pdb');
    const model = await plugin.builders.structure.createModel(trajectory);
    const structure = await plugin.builders.structure.createStructure(model);
    await plugin.builders.structure.representation.addRepresentation(structure, {
      type: 'spacefill',
      color: 'element-index',
    });
    const image = () => plugin.getImageRaw({ width: Size, height: Size });
    // pipelines compile in the background: frames drawn before they are ready skip draws
    await image();
    await plugin.canvas3d!.graphics.resources.pending.whenReady();
    return Array.from((await image()).data);
  } finally {
    plugin.dispose();
  }
}

// An item's theme layers together beyond the storage binding limit are read from up to 4 pages, each within it
describe.skipIf(!WebGPU.module)('data layers on WebGPU', () => {
  it('beyond the binding limit, read from several pages, draw the same image', async () => {
    const device = (WebGPU.module!.globals as { GPUDevice: { prototype: GPUDevice } }).GPUDevice.prototype;
    const createShaderModule = device.createShaderModule;
    const codes: string[] = [];
    device.createShaderModule = function (this: GPUDevice, descriptor: GPUShaderModuleDescriptor) {
      codes.push(descriptor.code);
      return createShaderModule.call(this, descriptor);
    };
    const warn = vi.spyOn(console, 'warn').mockImplementation(() => {});
    try {
      // the sizes of the item's layers, from the warning of a limit none fits
      await render(4);
      const sizes = /\(([\d, ]+) bytes\)/
        .exec(String(warn.mock.calls.find((c) => String(c[0]).includes('storage bindings'))?.[0]))![1]
        .split(', ')
        .map(Number);
      const largest = Math.max(...sizes);
      expect(sizes.reduce((a, b) => a + b, 0)).toBeGreaterThan(largest);
      warn.mockClear();

      const reference = await render();
      // not just the background
      expect(new Set(reference).size).toBeGreaterThan(2);
      // within the limit, one page: the programs of before
      expect(codes.some((c) => c.includes('sDataPage1'))).toBe(false);
      // the largest layer (positions) alone in the first page, the colors in the second
      codes.length = 0;
      expect(await render(largest)).toEqual(reference);
      expect(codes.some((c) => c.includes('sDataPage1'))).toBe(true);
      expect(warn).not.toHaveBeenCalled();

      // a layer alone beyond the limit: the item is not drawn, with a warning
      expect(await render(largest - 16)).not.toEqual(reference);
      expect(warn).toHaveBeenCalledWith(expect.stringContaining(`exceed 4 storage bindings of ${largest - 16} bytes`));
    } finally {
      device.createShaderModule = createShaderModule;
    }
  });
});
