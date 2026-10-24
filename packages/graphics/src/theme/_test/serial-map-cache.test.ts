/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Taylor Hoffmann <taylor@hoffmann.io>
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { beforeAll, describe, expect, it } from 'vitest';
import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import { CIF } from '@molstar/io/reader/cif';
import { trajectoryFromMmCIF } from '@molstar/model/formats/structure/mmcif';
import { Structure, StructureElement } from '@molstar/model/model/structure';
import type { ColorTheme } from '../color.js';
import { ChainIdColorThemeProvider } from '../color/chain-id.js';
import { EntityIdColorThemeProvider } from '../color/entity-id.js';
import { OperatorHklColorThemeProvider } from '../color/operator-hkl.js';

// two chains, A (entity 1) and B (entity 2), one atom each
const Cif = `data_TEST
loop_
_entity.id
_entity.type
1 polymer
2 polymer
loop_
_atom_site.group_PDB
_atom_site.id
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.label_entity_id
_atom_site.label_seq_id
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
_atom_site.auth_seq_id
_atom_site.auth_asym_id
_atom_site.pdbx_PDB_model_num
ATOM 1 C CA GLY A 1 1 0.0 0.0 0.0 1 A 1
ATOM 2 C CA ALA B 2 1 5.0 0.0 0.0 1 B 1
`;

function propertyKeys(structure: Structure) {
  return Object.keys(structure.currentPropertyData);
}

function createTheme<P extends PD.Params>(provider: ColorTheme.Provider<P, any>, structure: Structure) {
  const params = provider.getParams({ structure });
  return provider.factory({ structure }, PD.getDefaultValues(params));
}

function colors(theme: ColorTheme<any>, structure: Structure) {
  return structure.units.map((unit) =>
    theme.color(StructureElement.Location.create(structure, unit, unit.elements[0]), false),
  );
}

describe('color theme serial maps cached on the structure', () => {
  let structure: Structure;

  beforeAll(async () => {
    const parsed = await CIF.parseText(Cif).run();
    if (parsed.isError) throw new Error(parsed.message);
    const model = (await trajectoryFromMmCIF(parsed.result.blocks[0], parsed.result).run()).representative;
    structure = Structure.ofModel(model);
    expect(structure.units.length).toBe(2);
  });

  it('chain-id: a second theme reuses the map and colors the same', () => {
    const first = createTheme(ChainIdColorThemeProvider, structure);
    const keys = propertyKeys(structure).filter((k) => k.startsWith('__ChainIdColorTheme_'));
    expect(keys.length).toBe(1);
    const map = structure.currentPropertyData[keys[0]];

    const second = createTheme(ChainIdColorThemeProvider, structure);
    expect(structure.currentPropertyData[keys[0]]).toBe(map);
    expect(colors(second, structure)).toEqual(colors(first, structure));
    expect(new Set(colors(first, structure)).size).toBe(2);
  });

  it('chain-id: label and auth asym ids are cached apart', () => {
    const params = PD.getDefaultValues(ChainIdColorThemeProvider.getParams({ structure }));
    ChainIdColorThemeProvider.factory({ structure }, { ...params, asymId: 'label' });
    ChainIdColorThemeProvider.factory({ structure }, { ...params, asymId: 'auth' });
    const keys = propertyKeys(structure).filter((k) => k.startsWith('__ChainIdColorTheme_'));
    expect(keys.length).toBe(2);
  });

  it('entity-id: cached on the structure the theme is for, not on its root', () => {
    const sub = Structure.create([structure.units[1]], { parent: structure });
    const first = createTheme(EntityIdColorThemeProvider, sub);
    expect(propertyKeys(sub).some((k) => k.startsWith('__EntityIdColorTheme_'))).toBe(true);
    expect(propertyKeys(structure).some((k) => k.startsWith('__EntityIdColorTheme_'))).toBe(false);

    const second = createTheme(EntityIdColorThemeProvider, sub);
    expect(colors(second, sub)).toEqual(colors(first, sub));
  });

  it('operator-hkl: the params getter and the theme share one map', () => {
    OperatorHklColorThemeProvider.getParams({ structure });
    const keys = propertyKeys(structure).filter((k) => k.startsWith('__OperatorHklColorTheme_'));
    expect(keys.length).toBe(1);
    const data = structure.currentPropertyData[keys[0]];

    const theme = createTheme(OperatorHklColorThemeProvider, structure);
    expect(structure.currentPropertyData[keys[0]]).toBe(data);
    expect(colors(theme, structure).every((c) => !Number.isNaN(c))).toBe(true);
  });
});
