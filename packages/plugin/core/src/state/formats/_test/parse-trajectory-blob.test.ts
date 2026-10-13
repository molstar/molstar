/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import * as fs from 'node:fs';
import * as path from 'node:path';
import { PluginContext } from '@molstar/plugin/context';
import { PluginStateObject as SO, PluginStateTransform } from '@molstar/plugin/state/objects';
import { Mmcif, MmcifProvider, parseMmcifBlob } from '../trajectory/mmcif.js';
import { Pdb } from '../trajectory/pdb.js';

const mmcif = (id: string) =>
  [
    `data_${id}`,
    'loop_',
    '_atom_site.group_PDB',
    '_atom_site.id',
    '_atom_site.type_symbol',
    '_atom_site.label_atom_id',
    '_atom_site.label_comp_id',
    '_atom_site.label_asym_id',
    '_atom_site.label_entity_id',
    '_atom_site.label_seq_id',
    '_atom_site.Cartn_x',
    '_atom_site.Cartn_y',
    '_atom_site.Cartn_z',
    '_atom_site.pdbx_PDB_model_num',
    'ATOM 1 N N ALA A 1 1 0.0 0.0 0.0 1',
    '',
  ].join('\n');

const TestBlob = PluginStateTransform.BuiltIn({
  name: 'test-parse-trajectory-blob',
  display: 'Test blob',
  from: SO.Root,
  to: SO.Data.Blob,
})({
  apply: () =>
    new SO.Data.Blob([
      { id: '0', kind: 'string', data: mmcif('a') },
      { id: '1', kind: 'string', data: mmcif('b') },
    ]),
});

async function createPlugin(...formats: (typeof Mmcif)[]) {
  const plugin = new PluginContext({ behaviors: [], registry: [] });
  await plugin.init();
  plugin.dataFormats.clear();
  plugin.register(formats);
  const blob = await plugin.state.data.build().toRoot().apply(TestBlob).commit();
  return { plugin, blob };
}

const params = {
  formats: [
    { id: '0', format: 'cif' as const },
    { id: '1', format: 'cif' as const },
  ],
};

describe('parseTrajectory(blob)', () => {
  it('parses the blob through the registered mmCIF format', async () => {
    const { plugin, blob } = await createPlugin(Mmcif);
    expect(MmcifProvider.parseBlob).toBe(parseMmcifBlob);
    const trajectory = await plugin.builders.structure.parseTrajectory(blob, params);
    expect(trajectory.data!.frameCount).toBe(2);
  });

  it('fails with a clear error when mmCIF is not registered', async () => {
    const { plugin, blob } = await createPlugin(Pdb);
    await expect(plugin.builders.structure.parseTrajectory(blob, params)).rejects.toThrow(
      "parseTrajectory(blob) requires the 'mmcif' data format to be registered in this plugin.",
    );
  });

  it('is exported by the mmCIF module as a helper', async () => {
    const { plugin, blob } = await createPlugin();
    const trajectory = await parseMmcifBlob(plugin, blob, params);
    expect(trajectory.data!.frameCount).toBe(2);
  });

  it('the structure builder does not import the mmCIF module', () => {
    const source = fs.readFileSync(path.resolve(__dirname, '../../builder/structure.ts'), 'utf8');
    const imports = source.split('\n').filter((l) => l.startsWith('import ') && !l.startsWith('import type'));
    expect(imports.filter((l) => /formats\/(cif|trajectory)/.test(l))).toEqual([]);
  });
});
