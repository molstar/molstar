/**
 * Copyright (c) 2023-2024 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import type { Color } from '@molstar/core/util/color';
import type { Location } from '@molstar/model/model/location';
import type { ColorTheme } from '../color.js';
import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import type { ThemeDataContext } from '../theme.js';
import { ChainIdColorTheme, ChainIdColorThemeParams } from './chain-id.js';
import { UniformColorTheme, UniformColorThemeParams } from './uniform.js';
import { assertUnreachable } from '@molstar/core/util/type-helpers';
import { EntityIdColorTheme, EntityIdColorThemeParams } from './entity-id.js';
import { MoleculeTypeColorTheme, MoleculeTypeColorThemeParams } from './molecule-type.js';
import { EntitySourceColorTheme, EntitySourceColorThemeParams } from './entity-source.js';
import { ModelIndexColorTheme, ModelIndexColorThemeParams } from './model-index.js';
import { StructureIndexColorTheme, StructureIndexColorThemeParams } from './structure-index.js';
import { ColorThemeCategory } from './categories.js';
import { ResidueNameColorTheme, ResidueNameColorThemeParams } from './residue-name.js';
import type { ScaleLegend, TableLegend } from '@molstar/core/util/legend';
import { SecondaryStructureColorTheme, SecondaryStructureColorThemeParams } from './secondary-structure.js';
import { ElementSymbolColorTheme, ElementSymbolColorThemeParams } from './element-symbol.js';
import { TrajectoryIndexColorTheme, TrajectoryIndexColorThemeParams } from './trajectory-index.js';
import { hash2 } from '@molstar/core/data/util/hash-functions';
import { HydrophobicityColorTheme, HydrophobicityColorThemeParams } from './hydrophobicity.js';
import { UncertaintyColorTheme, UncertaintyColorThemeParams } from './uncertainty.js';
import { OccupancyColorTheme, OccupancyColorThemeParams } from './occupancy.js';
import { SequenceIdColorTheme, SequenceIdColorThemeParams } from './sequence-id.js';
import { PartialChargeColorTheme, PartialChargeColorThemeParams } from './partial-charge.js';

const Description = 'Uses separate themes for coloring mainchain and sidechain visuals.';

export const CartoonColorThemeParams = {
  mainchain: PD.MappedStatic('molecule-type', {
    uniform: PD.Group(UniformColorThemeParams),
    'chain-id': PD.Group(ChainIdColorThemeParams),
    'entity-id': PD.Group(EntityIdColorThemeParams),
    'entity-source': PD.Group(EntitySourceColorThemeParams),
    'molecule-type': PD.Group(MoleculeTypeColorThemeParams),
    'model-index': PD.Group(ModelIndexColorThemeParams),
    'structure-index': PD.Group(StructureIndexColorThemeParams),
    'secondary-structure': PD.Group(SecondaryStructureColorThemeParams),
    'trajectory-index': PD.Group(TrajectoryIndexColorThemeParams),
  }),
  sidechain: PD.MappedStatic('residue-name', {
    uniform: PD.Group(UniformColorThemeParams),
    'residue-name': PD.Group(ResidueNameColorThemeParams),
    'element-symbol': PD.Group(ElementSymbolColorThemeParams),
    hydrophobicity: PD.Group(HydrophobicityColorThemeParams),
    uncertainty: PD.Group(UncertaintyColorThemeParams),
    occupancy: PD.Group(OccupancyColorThemeParams),
    'sequence-id': PD.Group(SequenceIdColorThemeParams),
    'partial-charge': PD.Group(PartialChargeColorThemeParams),
  }),
};
export type CartoonColorThemeParams = typeof CartoonColorThemeParams;
export function getCartoonColorThemeParams(ctx: ThemeDataContext) {
  const params = PD.clone(CartoonColorThemeParams);
  return params;
}

type CartoonColorThemeProps = PD.Values<CartoonColorThemeParams>;

function getMainchainTheme(ctx: ThemeDataContext, props: CartoonColorThemeProps['mainchain']) {
  switch (props.name) {
    case 'uniform':
      return UniformColorTheme(ctx, props.params);
    case 'chain-id':
      return ChainIdColorTheme(ctx, props.params);
    case 'entity-id':
      return EntityIdColorTheme(ctx, props.params);
    case 'entity-source':
      return EntitySourceColorTheme(ctx, props.params);
    case 'molecule-type':
      return MoleculeTypeColorTheme(ctx, props.params);
    case 'model-index':
      return ModelIndexColorTheme(ctx, props.params);
    case 'structure-index':
      return StructureIndexColorTheme(ctx, props.params);
    case 'secondary-structure':
      return SecondaryStructureColorTheme(ctx, props.params);
    case 'trajectory-index':
      return TrajectoryIndexColorTheme(ctx, props.params);
    default:
      assertUnreachable(props);
  }
}

function getSidechainTheme(ctx: ThemeDataContext, props: CartoonColorThemeProps['sidechain']) {
  switch (props.name) {
    case 'uniform':
      return UniformColorTheme(ctx, props.params);
    case 'residue-name':
      return ResidueNameColorTheme(ctx, props.params);
    case 'element-symbol':
      return ElementSymbolColorTheme(ctx, props.params);
    case 'hydrophobicity':
      return HydrophobicityColorTheme(ctx, props.params);
    case 'uncertainty':
      return UncertaintyColorTheme(ctx, props.params);
    case 'occupancy':
      return OccupancyColorTheme(ctx, props.params);
    case 'sequence-id':
      return SequenceIdColorTheme(ctx, props.params);
    case 'partial-charge':
      return PartialChargeColorTheme(ctx, props.params);
    default:
      assertUnreachable(props);
  }
}

export function CartoonColorTheme(
  ctx: ThemeDataContext,
  props: PD.Values<CartoonColorThemeParams>,
): ColorTheme<CartoonColorThemeParams> {
  const mainchain = getMainchainTheme(ctx, props.mainchain);
  const sidechain = getSidechainTheme(ctx, props.sidechain);

  const contextHash = hash2(mainchain.contextHash ?? 0, sidechain.contextHash ?? 0);

  function color(location: Location, isSecondary: boolean): Color {
    return isSecondary ? mainchain.color(location, false) : sidechain.color(location, false);
  }

  let legend: ScaleLegend | TableLegend | undefined = mainchain.legend;
  if (mainchain.legend?.kind === 'table-legend' && sidechain.legend?.kind === 'table-legend') {
    legend = {
      kind: 'table-legend',
      table: [...mainchain.legend.table, ...sidechain.legend.table],
    };
  }

  return {
    factory: CartoonColorTheme,
    granularity: 'group',
    preferSmoothing: false,
    color,
    props,
    contextHash,
    description: Description,
    legend,
  };
}

export const CartoonColorThemeProvider: ColorTheme.Provider<CartoonColorThemeParams, 'cartoon'> = {
  name: 'cartoon',
  label: 'Cartoon',
  category: ColorThemeCategory.Misc,
  factory: CartoonColorTheme,
  getParams: getCartoonColorThemeParams,
  defaultValues: PD.getDefaultValues(CartoonColorThemeParams),
  isApplicable: (ctx: ThemeDataContext) => !!ctx.structure,
};
