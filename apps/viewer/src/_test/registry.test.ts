/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import * as fs from 'node:fs';
import * as path from 'node:path';
import { PluginContext } from '@molstar/plugin/context';
import { DefaultFormats, DefaultRegistry } from '@molstar/plugin/default-registry';
import type { DataFormatProvider } from '@molstar/plugin/state/formats/provider';
import { DefaultPluginSpec } from '@molstar/plugin/default-spec';
import { createViewerRegistry } from '@molstar/viewer/registry';

const mmcif = DefaultFormats.formats!.find((f) => f.name === 'mmcif')!;

function customProvider(label: string): DataFormatProvider.Unnamed {
  return { ...mmcif, label, name: undefined };
}

async function createPlugin(customFormats?: [string, DataFormatProvider.Unnamed][]) {
  const registry = createViewerRegistry(customFormats);
  const plugin = new PluginContext({ ...DefaultPluginSpec(), registry });
  await plugin.init();
  return { registry, plugin };
}

describe('createViewerRegistry', () => {
  it('lists DefaultRegistry followed by the custom formats entry', async () => {
    const { registry, plugin } = await createPlugin([['extra', customProvider('Custom extra')]]);
    expect(registry.slice(0, -1)).toEqual(DefaultRegistry);
    expect(registry[registry.length - 1].formats!.map((f) => f.name)).toEqual(['extra']);
    expect(plugin.managers.animation.animations[0].name).toBe('built-in.animate-model-index');

    // the custom provider is registered after every built-in one
    expect(plugin.dataFormats.list[plugin.dataFormats.list.length - 1].name).toBe('extra');
    plugin.dispose();
  });

  it('registers a custom format under its name, which may be new', async () => {
    const { plugin } = await createPlugin([['extra', customProvider('Custom extra')]]);
    expect(plugin.dataFormats.get('extra')).toMatchObject({ name: 'extra', label: 'Custom extra' });
    expect(plugin.dataFormats.has('mmcif')).toBe(true);
    plugin.dispose();
  });

  it('lets a custom format override a built-in one, as in 5.x', async () => {
    const { registry, plugin } = await createPlugin([['mmcif', customProvider('Custom mmcif')]]);

    const names = plugin.dataFormats.list.map((f) => f.name);
    expect(names.filter((n) => n === 'mmcif')).toHaveLength(1);
    expect(plugin.dataFormats.get('mmcif')).toMatchObject({ name: 'mmcif', label: 'Custom mmcif' });
    expect(plugin.dataFormats.get('mmcif')).not.toBe(mmcif);

    // the default entry is replaced by a copy without the overridden provider, and is not mutated
    expect(DefaultFormats.formats!.includes(mmcif)).toBe(true);
    const formatEntries = registry.filter((e) => e.formats);
    expect(formatEntries.some((e) => e === DefaultFormats)).toBe(false);
    expect(formatEntries[0].formats!.some((f) => f.name === 'mmcif')).toBe(false);
    plugin.dispose();
  });

  it('loads data with the custom provider of an overridden built-in format', async () => {
    const pdb = DefaultFormats.formats!.find((f) => f.name === 'pdb')!;
    const parse = jest.fn((...args: Parameters<typeof pdb.parse>) => pdb.parse(...args));
    const custom = { ...pdb, label: 'Custom pdb', name: undefined, parse };
    const { plugin } = await createPlugin([['pdb', custom]]);
    expect(plugin.dataFormats.get('pdb').label).toBe('Custom pdb');

    const text = fs.readFileSync(path.resolve(__dirname, '../../../../smoke/fixtures/tiny.pdb'), 'utf8');
    const data = await plugin.builders.data.rawData({ data: text });
    const trajectory = await plugin.builders.structure.parseTrajectory(data, 'pdb');
    expect(parse).toHaveBeenCalledTimes(1);
    expect(trajectory.data!.frameCount).toBe(1);

    // the provider is also found by the format auto-detection of the data
    expect(plugin.dataFormats.auto({ name: 'tiny.pdb', ext: 'pdb' } as any, data.cell!.obj as any)).toBe(
      plugin.dataFormats.get('pdb'),
    );
    plugin.dispose();
  });
});
