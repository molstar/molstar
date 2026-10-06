/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import * as fs from 'node:fs';
import * as path from 'node:path';
import { Script } from '@molstar/model/script/script';
import { StructureElement } from '@molstar/model/model/structure';
import { PluginContext } from '@molstar/plugin/context';
import { DefaultPluginSpec } from '@molstar/plugin/default-spec';

const crambin = fs.readFileSync(path.resolve(__dirname, '../../../../../data/examples/1crn.cif'), 'utf8');

describe('the default plugin spec and script languages', () => {
  it('evaluates PyMOL, VMD and Jmol scripts (the default spec imports transpilers/all)', async () => {
    const plugin = new PluginContext(DefaultPluginSpec());
    await plugin.init();
    const data = await plugin.builders.data.rawData({ data: crambin });
    const trajectory = await plugin.builders.structure.parseTrajectory(data, 'mmcif');
    const model = await plugin.builders.structure.createModel(trajectory);
    const structureRef = await plugin.builders.structure.createStructure(model);
    const structure = structureRef.obj!.data;

    expect(Script.getAvailableLanguages().sort()).toEqual(['jmol', 'mol-script', 'pymol', 'vmd']);

    const pymol = Script.toLoci({ language: 'pymol', expression: 'resn ALA' }, structure);
    const vmd = Script.toLoci({ language: 'vmd', expression: 'resname ALA' }, structure);
    const molScript = Script.toLoci(
      { language: 'mol-script', expression: '(sel.atom.atom-groups :residue-test (= atom.label_comp_id ALA))' },
      structure,
    );
    expect(StructureElement.Loci.size(pymol)).toBeGreaterThan(0);
    expect(StructureElement.Loci.size(pymol)).toBe(StructureElement.Loci.size(molScript));
    expect(StructureElement.Loci.size(vmd)).toBe(StructureElement.Loci.size(pymol));
    plugin.dispose();
  });
});
