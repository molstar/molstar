/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import * as fs from 'node:fs';
import * as path from 'node:path';
import type { PluginContext } from '@molstar/plugin/context';
import { nucleic, protein } from '@molstar/plugin/state/queries/structure/type';
import { presetSelectionComponent } from '../representation-presets/types.js';

describe('presetSelectionComponent', () => {
  it('creates the component from the given query under the tag', async () => {
    const tryCreateComponentFromSelection = jest.fn(async () => undefined);
    const plugin = { builders: { structure: { tryCreateComponentFromSelection } } } as unknown as PluginContext;
    const params = { label: 'Protein', tags: ['x'] };
    const structure = {} as any;

    await presetSelectionComponent(plugin, structure, protein, 'protein', params);
    await presetSelectionComponent(plugin, structure, nucleic, 'nucleic');
    expect(tryCreateComponentFromSelection.mock.calls).toEqual([
      [structure, protein, 'selection-protein', params],
      [structure, nucleic, 'selection-nucleic', undefined],
    ]);
  });

  it('does not import the selection query catalog', () => {
    const source = fs.readFileSync(path.resolve(__dirname, '../representation-presets/types.ts'), 'utf8');
    const imports = source.split('\n').filter((l) => l.startsWith('import ') && !l.startsWith('import type'));
    expect(imports.filter((l) => /queries\/structure\/catalog/.test(l))).toEqual([]);
  });
});
