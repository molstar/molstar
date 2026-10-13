/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import * as fs from 'node:fs';
import * as path from 'node:path';
import { PluginContext } from '@molstar/plugin/context';
import { DefaultPluginSpec } from '@molstar/plugin/default-spec';

const crambin = fs.readFileSync(path.resolve(__dirname, '../../../../../../../../data/examples/1crn.cif'), 'utf8');

async function createPluginWithStructure() {
  const plugin = new PluginContext(DefaultPluginSpec());
  await plugin.init();
  const data = await plugin.builders.data.rawData({ data: crambin });
  const trajectory = await plugin.builders.structure.parseTrajectory(data, 'mmcif');
  const model = await plugin.builders.structure.createModel(trajectory);
  const structure = await plugin.builders.structure.createStructure(model);
  return { plugin, structure };
}

describe('building a representation with unregistered names in a plugin with the default spec', () => {
  it('warns naming the provider and scope and uses the registry default', async () => {
    const { plugin, structure } = await createPluginWithStructure();
    const registry = plugin.representation.structure.registry;
    const warn = jest.spyOn(plugin.log, 'warn').mockImplementation(() => {});

    const repr = await plugin.builders.structure.representation.addRepresentation(structure, {
      type: 'no-such-repr' as any,
      color: 'no-such-color' as any,
      size: 'no-such-size' as any,
    });
    const messages = warn.mock.calls.map((c) => String(c[0]));
    expect(messages).toEqual(
      expect.arrayContaining([
        "Structure representation 'no-such-repr' is not registered in this plugin; the registry default is used",
        "Structure color theme 'no-such-color' is not registered in this plugin; the registry default is used",
        "Structure size theme 'no-such-size' is not registered in this plugin; the registry default is used",
      ]),
    );

    const params = repr.cell!.transform.params as any;
    const defaultRepr = registry.default!.provider;
    expect(params.type.name).toBe(registry.default!.name);
    expect(params.colorTheme.name).toBe(defaultRepr.defaultColorTheme.name);
    expect(params.sizeTheme.name).toBe(defaultRepr.defaultSizeTheme.name);
    expect(repr.cell!.status).toBe('ok');
    plugin.dispose();
  });

  it('does not warn for registered names', async () => {
    const { plugin, structure } = await createPluginWithStructure();
    const warn = jest.spyOn(plugin.log, 'warn').mockImplementation(() => {});
    const repr = await plugin.builders.structure.representation.addRepresentation(structure, {
      type: 'ball-and-stick',
      color: 'chain-id',
    });
    expect(warn).not.toHaveBeenCalled();
    expect((repr.cell!.transform.params as any).type.name).toBe('ball-and-stick');
    plugin.dispose();
  });
});
