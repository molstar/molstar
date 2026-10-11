/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { HeadlessPluginContext } from '@molstar/plugin-headless/context';
import type { ExternalModules } from '@molstar/plugin-headless/screenshot';
import type { WebGPUContext } from '@molstar/graphics/gl/webgpu/context';
import { DefaultPluginSpec } from '@molstar/plugin/spec';

/**
 * The optional `webgpu` package (Dawn's Node.js bindings), with a usable adapter. It is an ES module, so Jest needs
 * `--experimental-vm-modules` to load it.
 */
async function getWebGPUModule(): Promise<ExternalModules['webgpu'] | undefined> {
    let module: ExternalModules['webgpu'];
    try {
        module = (await import('webgpu' as string)) as ExternalModules['webgpu'];
    } catch (e) {
        console.warn(`Skipping headless WebGPU test: optional package "webgpu" is unavailable (${(e as Error).message}).`);
        return;
    }
    const adapter = await module!.create([]).requestAdapter();
    if (!adapter) {
        console.warn('Skipping headless WebGPU test: no WebGPU adapter.');
        return;
    }
    return module;
}

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

describe('HeadlessPluginContext on WebGPU', () => {
    it('renders a structure', async () => {
        const webgpu = await getWebGPUModule();
        if (!webgpu) return;

        const plugin = await HeadlessPluginContext.create(
            { webgpu },
            DefaultPluginSpec(),
            { width: 64, height: 64 },
            { backend: 'webgpu' },
        );
        try {
            await plugin.init();
            expect(plugin.canvas3d!.graphics.capabilities.backend).toBe('webgpu');

            const data = await plugin.builders.data.rawData({ data: Pdb, label: 'test' });
            const trajectory = await plugin.builders.structure.parseTrajectory(data, 'pdb');
            await plugin.builders.structure.hierarchy.applyPreset(trajectory, 'default');

            const image = await plugin.getImageRaw({ width: 64, height: 64 });
            expect(image.data.length).toBe(64 * 64 * 4);
            // the structure covers part of the image, the rest is background
            const first = image.data.slice(0, 4).join();
            let differing = 0;
            for (let i = 0; i < image.data.length; i += 4) {
                if (image.data.slice(i, i + 4).join() !== first) differing += 1;
            }
            expect(differing).toBeGreaterThan(64);
            expect(differing).toBeLessThan(64 * 64);

            // device loss: the context replaces the device and recreates its resources, the image is the same
            const graphics = plugin.canvas3d!.graphics as WebGPUContext;
            const lostDevice = graphics.device;
            const restored = new Promise<void>((resolve) => graphics.contextRestored.subscribe(() => resolve()));
            lostDevice.destroy();
            await restored;
            expect(graphics.device).not.toBe(lostDevice);
            expect(graphics.isContextLost).toBe(false);
            const restoredImage = await plugin.getImageRaw({ width: 64, height: 64 });
            expect(restoredImage.data).toEqual(image.data);
        } finally {
            plugin.dispose();
        }
    });
});
