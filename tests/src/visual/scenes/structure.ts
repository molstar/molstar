/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Column, Table } from '@molstar/core/data/db';
import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import { RuntimeContext } from '@molstar/core/task';
import { type GraphicsContext } from '@molstar/graphics/gl/context';
import { Representation, type RepresentationProvider } from '@molstar/graphics/repr/representation';
import { ColorTheme } from '@molstar/graphics/theme/color';
import { BuiltInColorThemes } from '@molstar/graphics/theme/color/catalog';
import { SizeTheme } from '@molstar/graphics/theme/size';
import { BuiltInSizeThemes } from '@molstar/graphics/theme/size/catalog';
import { Theme } from '@molstar/graphics/theme/theme';
import { createModels } from '@molstar/model/formats/structure/basic/parser';
import { BasicSchema, createBasic } from '@molstar/model/formats/structure/basic/schema';
import { ComponentBuilder } from '@molstar/model/formats/structure/common/component';
import { EntityBuilder } from '@molstar/model/formats/structure/common/entity';
import { Structure } from '@molstar/model/model/structure';
import { CIF } from '@molstar/io/reader/cif';
import { trajectoryFromMmCIF } from '@molstar/model/formats/structure/mmcif';
import { MoleculeType } from '@molstar/model/model/structure/model/types';

let themeRegistries: { colorThemeRegistry: ColorTheme.Registry; sizeThemeRegistry: SizeTheme.Registry } | undefined;
/** The built-in color and size themes */
export function getThemeRegistries() {
  if (!themeRegistries) {
    // registries start empty
    const colorThemeRegistry = ColorTheme.createRegistry();
    for (const p of Object.values(BuiltInColorThemes)) colorThemeRegistry.add(p as ColorTheme.Provider<any, any, any>);
    const sizeThemeRegistry = SizeTheme.createRegistry();
    for (const p of Object.values(BuiltInSizeThemes)) sizeThemeRegistry.add(p as SizeTheme.Provider<any, any>);
    themeRegistries = { colorThemeRegistry, sizeThemeRegistry };
  }
  return themeRegistries;
}

/** Two chains of carbon atoms on helices around the z axis, one atom per residue */
export async function helixStructure(atomsPerChain: number): Promise<Structure> {
  const chains = ['A', 'B'];
  const count = chains.length * atomsPerChain;
  const x = new Float32Array(count),
    y = new Float32Array(count),
    z = new Float32Array(count);
  const chainOf = new Int32Array(count);
  for (let c = 0; c < chains.length; ++c) {
    for (let i = 0; i < atomsPerChain; ++i) {
      const k = c * atomsPerChain + i;
      const a = i * 0.6 + c * Math.PI;
      x[k] = 4 * Math.cos(a);
      y[k] = 4 * Math.sin(a);
      z[k] = i * 0.5 - atomsPerChain * 0.25;
      chainOf[k] = c;
    }
  }
  return atomStructure('helices', x, y, z, chains, chainOf);
}

/**
 * A structure of carbon atoms ('CA' of 'ALA' residues, one atom per residue) at the given coordinates, atom `i` in
 * chain `chains[chainOf[i]]`; the atoms of a chain must be consecutive
 */
export async function atomStructure(
  name: string,
  x: Float32Array,
  y: Float32Array,
  z: Float32Array,
  chains: readonly string[],
  chainOf: Int32Array,
): Promise<Structure> {
  const count = x.length;
  const asymIds: string[] = [];
  const seqIds: number[] = [];
  let seqId = 0;
  for (let k = 0; k < count; ++k) {
    if (k === 0 || chainOf[k] !== chainOf[k - 1]) seqId = 0;
    asymIds.push(chains[chainOf[k]]);
    seqIds.push(++seqId);
  }

  const auth_asym_id = Column.ofStringArray(asymIds);
  const auth_atom_id = Column.ofConst('CA', count, Column.Schema.str);
  const auth_comp_id = Column.ofConst('ALA', count, Column.Schema.str);
  const auth_seq_id = Column.ofIntArray(seqIds);
  const atom_site = Table.ofPartialColumns(
    BasicSchema.atom_site,
    {
      auth_asym_id,
      auth_atom_id,
      auth_comp_id,
      auth_seq_id,
      Cartn_x: Column.ofFloatArray(x),
      Cartn_y: Column.ofFloatArray(y),
      Cartn_z: Column.ofFloatArray(z),
      id: Column.ofIntArray(seqIds.map((_, i) => i + 1)),
      label_asym_id: auth_asym_id,
      label_atom_id: auth_atom_id,
      label_comp_id: auth_comp_id,
      label_seq_id: auth_seq_id,
      label_entity_id: Column.ofConst('1', count, Column.Schema.str),
      occupancy: Column.ofConst(1, count, Column.Schema.float),
      type_symbol: Column.ofConst('C', count, Column.Schema.str),
      pdbx_PDB_model_num: Column.ofConst(1, count, Column.Schema.int),
    },
    count,
  );

  const entityBuilder = new EntityBuilder();
  for (const chain of chains) entityBuilder.getEntityId('ALA', MoleculeType.Protein, chain);
  const componentBuilder = new ComponentBuilder(auth_seq_id, auth_atom_id);
  componentBuilder.add('ALA', 0);

  const basic = createBasic({
    entity: entityBuilder.getEntityTable(),
    chem_comp: componentBuilder.getChemCompTable(),
    atom_site,
  });
  const trajectory = await createModels(
    basic,
    { kind: 'visual-test', name, data: undefined },
    RuntimeContext.Synchronous,
  );
  return Structure.ofModel(trajectory.representative);
}

export async function structureRepr<P extends PD.Params>(
  graphics: GraphicsContext,
  provider: RepresentationProvider<Structure, P>,
  structure: Structure,
  props: Partial<PD.Values<P>>,
  colorTheme?: { name: string; params?: {} },
): Promise<Representation.Any> {
  const registries = getThemeRegistries();
  const repr = provider.factory({ graphics, ...registries }, provider.getParams);
  repr.setTheme(
    Theme.create(
      registries,
      { structure },
      {
        colorTheme: {
          name: colorTheme?.name ?? provider.defaultColorTheme.name,
          params: colorTheme?.params ?? provider.defaultColorTheme.props ?? {},
        },
        sizeTheme: { name: provider.defaultSizeTheme.name, params: provider.defaultSizeTheme.props ?? {} },
      },
    ),
  );
  await repr.createOrUpdate(props, structure).run();
  return repr as unknown as Representation.Any;
}

/**
 * A PDB entry of `data/examples` (served with the repository), its first model: units per chain, as representations of
 * entries have them (unlike the globules, chains have cavities inside, adding to their surfaces)
 */
export async function exampleStructure(file: string): Promise<Structure> {
  const response = await fetch(`/data/examples/${file}`);
  if (!response.ok) throw new Error(`failed to load '${file}': ${response.status}`);
  const parsed = await CIF.parse(await response.text()).run();
  if (parsed.isError) throw new Error(`failed to parse '${file}': ${parsed.message}`);
  const trajectory = await trajectoryFromMmCIF(parsed.result.blocks[0]).run();
  return Structure.ofModel(trajectory.representative);
}
