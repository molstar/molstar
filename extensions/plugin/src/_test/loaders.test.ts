/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import * as fs from 'node:fs';
import * as path from 'node:path';
import { PluginContext } from '@molstar/plugin/context';
import { DefaultPluginSpec } from '@molstar/plugin/default-spec';
import { loadTrajectory } from '@molstar/plugin-extension/loaders';

const pdb = fs.readFileSync(path.resolve(__dirname, '../../../../smoke/fixtures/tiny.pdb'), 'utf8');

function lammpstrj(frames: number) {
  const out: string[] = [];
  for (let f = 0; f < frames; f++) {
    out.push('ITEM: TIMESTEP', `${f}`, 'ITEM: NUMBER OF ATOMS', '3', 'ITEM: BOX BOUNDS pp pp pp');
    out.push('0 10', '0 10', '0 10', 'ITEM: ATOMS id type x y z');
    for (let i = 0; i < 3; i++) out.push(`${i + 1} 1 ${i * 1.45 + f} 0 0`);
  }
  return out.join('\n') + '\n';
}

const modelAndCoords = (preset?: 'default' | 'all-models') => ({
  model: { kind: 'model-data' as const, data: pdb, format: 'pdb' as const },
  coordinates: { kind: 'coordinates-data' as const, data: lammpstrj(2), format: 'lammpstrj' as const },
  preset,
});

async function createPlugin() {
  const plugin = new PluginContext(DefaultPluginSpec());
  await plugin.init();
  return plugin;
}

describe('loadTrajectory with the default plugin', () => {
  it("applies the 'all-models' preset through its alias, with one structure per frame", async () => {
    const plugin = await createPlugin();
    expect(plugin.builders.structure.hierarchy.has('all-models')).toBe(true);
    const { preset } = await loadTrajectory(plugin, modelAndCoords('all-models'));
    expect(preset.structures).toHaveLength(2);
    expect(preset.models).toHaveLength(2);
    plugin.dispose();
  });

  it("applies the 'default' preset by default", async () => {
    const plugin = await createPlugin();
    const { preset } = await loadTrajectory(plugin, modelAndCoords());
    expect(preset.structure).toBeTruthy();
    expect(preset.structures).toBeUndefined();
    plugin.dispose();
  });

  it('fails with a clear error for an unknown preset', async () => {
    const plugin = await createPlugin();
    await expect(loadTrajectory(plugin, modelAndCoords('no-such-preset' as any))).rejects.toThrow(/no-such-preset/);
    plugin.dispose();
  });
});

describe('loadTrajectory without the format', () => {
  it('fails with a clear error naming the unregistered format', async () => {
    const plugin = new PluginContext({ ...DefaultPluginSpec(), registry: [] });
    await plugin.init();
    await expect(loadTrajectory(plugin, modelAndCoords())).rejects.toThrow(/'pdb' is not a supported data format/);
    plugin.dispose();
  });
});
