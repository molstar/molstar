/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { PluginContext } from '@molstar/plugin/context';
import type { PluginRegistryEntry } from '@molstar/plugin/spec';
import { Ccp4, Ccp4Provider } from '../volume/ccp4.js';

/** A little-endian CCP4 map of a gaussian blob, `n` grid points along each axis of a cubic cell. */
function createCcp4(n: number) {
  const buffer = new ArrayBuffer(1024 + n * n * n * 4);
  const view = new DataView(buffer);
  const int = (word: number, value: number) => view.setInt32(word * 4, value, true);
  const float = (word: number, value: number) => view.setFloat32(word * 4, value, true);

  const values = new Float32Array(buffer, 1024, n * n * n);
  let min = Infinity;
  let max = -Infinity;
  let sum = 0;
  for (let z = 0; z < n; z++) {
    for (let y = 0; y < n; y++) {
      for (let x = 0; x < n; x++) {
        const d = (x - n / 2) ** 2 + (y - n / 2) ** 2 + (z - n / 2) ** 2;
        const v = Math.exp(-d / 4);
        values[(z * n + y) * n + x] = v;
        min = Math.min(min, v);
        max = Math.max(max, v);
        sum += v;
      }
    }
  }

  int(0, n); // NC
  int(1, n); // NR
  int(2, n); // NS
  int(3, 2); // MODE: float32
  int(7, n); // NX
  int(8, n); // NY
  int(9, n); // NZ
  for (let i = 0; i < 3; i++) float(10 + i, n); // cell lengths
  for (let i = 0; i < 3; i++) float(13 + i, 90); // cell angles
  int(16, 1); // MAPC
  int(17, 2); // MAPR
  int(18, 3); // MAPS
  float(19, min);
  float(20, max);
  float(21, sum / values.length);
  int(22, 1); // ISPG
  // 'MAP ' and the little-endian machine stamp
  [...'MAP '].forEach((c, i) => view.setUint8(52 * 4 + i, c.charCodeAt(0)));
  view.setUint8(53 * 4, 68);
  view.setUint8(53 * 4 + 1, 65);
  return new Uint8Array(buffer);
}

/** A plugin whose volume representation and theme registries are empty, so that only the entries provide anything. */
async function createPlugin(registry: PluginRegistryEntry[]) {
  const plugin = new PluginContext({ behaviors: [], registry });
  const { volume } = plugin.representation;
  volume.registry.clear();
  volume.themes.colorThemeRegistry.clear();
  volume.themes.sizeThemeRegistry.clear();

  const warnings: string[] = [];
  plugin.events.log.subscribe((e) => {
    if (e.type === 'warning') warnings.push(e.message);
  });
  await plugin.init();
  return { plugin, warnings };
}

describe('loading a volume with only its format entry registered', () => {
  it('shows the isosurface without unregistered-name warnings', async () => {
    const { plugin, warnings } = await createPlugin([Ccp4]);
    expect(plugin.representation.volume.registry.has('isosurface')).toBe(true);
    expect(plugin.representation.volume.registry.has('direct-volume')).toBe(false);
    expect(plugin.representation.volume.themes.colorThemeRegistry.has('uniform')).toBe(true);
    expect(plugin.representation.volume.themes.sizeThemeRegistry.has('uniform')).toBe(true);

    const data = await plugin.builders.data.rawData(
      { data: createCcp4(8), label: 'blob.ccp4' },
      { state: { isGhost: true } },
    );
    const parsed = await Ccp4Provider.parse(plugin, data);
    const visuals = await Ccp4Provider.visuals!(plugin, parsed);

    expect(visuals.length).toBe(1);
    const repr = visuals[0].cell?.obj;
    expect(repr?.type.name).toBe('Volume 3D');
    expect(visuals[0].cell?.transform.params).toMatchObject({
      type: { name: 'isosurface' },
      colorTheme: { name: 'uniform' },
      sizeTheme: { name: 'uniform' },
    });
    expect(warnings).toEqual([]);
    await plugin.dispose();
  });

  it('fails to show anything when the entry is not registered', async () => {
    const { plugin } = await createPlugin([]);
    const data = await plugin.builders.data.rawData(
      { data: createCcp4(8), label: 'blob.ccp4' },
      { state: { isGhost: true } },
    );
    const parsed = await Ccp4Provider.parse(plugin, data);
    await expect(Ccp4Provider.visuals!(plugin, parsed)).rejects.toThrow();
    await plugin.dispose();
  });
});
