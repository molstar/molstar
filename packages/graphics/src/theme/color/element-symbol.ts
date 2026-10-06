/**
 * Copyright (c) 2018-2024 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author Adam Midlik <midlik@gmail.com>
 * @author Hakan Akgül <hakan-akgul@outlook.com>
 */

import type { ElementSymbol } from '@molstar/model/model/structure/model/types';
import { Color, ColorMap } from '@molstar/core/util/color';
import { StructureElement, Unit, Bond } from '@molstar/model/model/structure';
import type { Location } from '@molstar/model/model/location';
import type { ColorTheme } from '../color.js';
import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import type { ThemeDataContext } from '../theme.js';
import { TableLegend } from '@molstar/core/util/legend';
import { getAdjustedColorMap } from '@molstar/core/util/color/color';
import { getColorMapParams } from '@molstar/core/util/color/params';
import { ChainIdColorTheme, ChainIdColorThemeParams } from './chain-id.js';
import { OperatorNameColorThemeParams, OperatorNameColorTheme } from './operator-name.js';
import { EntityIdColorTheme, EntityIdColorThemeParams } from './entity-id.js';
import { assertUnreachable } from '@molstar/core/util/type-helpers';
import { EntitySourceColorTheme, EntitySourceColorThemeParams } from './entity-source.js';
import { ModelIndexColorTheme, ModelIndexColorThemeParams } from './model-index.js';
import { StructureIndexColorTheme, StructureIndexColorThemeParams } from './structure-index.js';
import { ColorThemeCategory } from './categories.js';
import { UnitIndexColorTheme, UnitIndexColorThemeParams } from './unit-index.js';
import { UniformColorTheme, UniformColorThemeParams } from './uniform.js';
import { TrajectoryIndexColorTheme, TrajectoryIndexColorThemeParams } from './trajectory-index.js';

// from Jmol http://jmol.sourceforge.net/jscolors/ (or 0xFFFFFF)
export const ElementSymbolColors = ColorMap({
  H: 0xffffff,
  D: 0xffffc0,
  T: 0xffffa0,
  HE: 0xd9ffff,
  LI: 0xcc80ff,
  BE: 0xc2ff00,
  B: 0xffb5b5,
  C: 0x909090,
  N: 0x3050f8,
  O: 0xff0d0d,
  F: 0x90e050,
  NE: 0xb3e3f5,
  NA: 0xab5cf2,
  MG: 0x8aff00,
  AL: 0xbfa6a6,
  SI: 0xf0c8a0,
  P: 0xff8000,
  S: 0xffff30,
  CL: 0x1ff01f,
  AR: 0x80d1e3,
  K: 0x8f40d4,
  CA: 0x3dff00,
  SC: 0xe6e6e6,
  TI: 0xbfc2c7,
  V: 0xa6a6ab,
  CR: 0x8a99c7,
  MN: 0x9c7ac7,
  FE: 0xe06633,
  CO: 0xf090a0,
  NI: 0x50d050,
  CU: 0xc88033,
  ZN: 0x7d80b0,
  GA: 0xc28f8f,
  GE: 0x668f8f,
  AS: 0xbd80e3,
  SE: 0xffa100,
  BR: 0xa62929,
  KR: 0x5cb8d1,
  RB: 0x702eb0,
  SR: 0x00ff00,
  Y: 0x94ffff,
  ZR: 0x94e0e0,
  NB: 0x73c2c9,
  MO: 0x54b5b5,
  TC: 0x3b9e9e,
  RU: 0x248f8f,
  RH: 0x0a7d8c,
  PD: 0x006985,
  AG: 0xc0c0c0,
  CD: 0xffd98f,
  IN: 0xa67573,
  SN: 0x668080,
  SB: 0x9e63b5,
  TE: 0xd47a00,
  I: 0x940094,
  XE: 0x940094,
  CS: 0x57178f,
  BA: 0x00c900,
  LA: 0x70d4ff,
  CE: 0xffffc7,
  PR: 0xd9ffc7,
  ND: 0xc7ffc7,
  PM: 0xa3ffc7,
  SM: 0x8fffc7,
  EU: 0x61ffc7,
  GD: 0x45ffc7,
  TB: 0x30ffc7,
  DY: 0x1fffc7,
  HO: 0x00ff9c,
  ER: 0x00e675,
  TM: 0x00d452,
  YB: 0x00bf38,
  LU: 0x00ab24,
  HF: 0x4dc2ff,
  TA: 0x4da6ff,
  W: 0x2194d6,
  RE: 0x267dab,
  OS: 0x266696,
  IR: 0x175487,
  PT: 0xd0d0e0,
  AU: 0xffd123,
  HG: 0xb8b8d0,
  TL: 0xa6544d,
  PB: 0x575961,
  BI: 0x9e4fb5,
  PO: 0xab5c00,
  AT: 0x754f45,
  RN: 0x428296,
  FR: 0x420066,
  RA: 0x007d00,
  AC: 0x70abfa,
  TH: 0x00baff,
  PA: 0x00a1ff,
  U: 0x008fff,
  NP: 0x0080ff,
  PU: 0x006bff,
  AM: 0x545cf2,
  CM: 0x785ce3,
  BK: 0x8a4fe3,
  CF: 0xa136d4,
  ES: 0xb31fd4,
  FM: 0xb31fba,
  MD: 0xb30da6,
  NO: 0xbd0d87,
  LR: 0xc70066,
  RF: 0xcc0059,
  DB: 0xd1004f,
  SG: 0xd90045,
  BH: 0xe00038,
  HS: 0xe6002e,
  MT: 0xeb0026,
  DS: 0xffffff,
  RG: 0xffffff,
  CN: 0xffffff,
  UUT: 0xffffff,
  NH: 0xffffff,
  FL: 0xffffff,
  UUP: 0xffffff,
  MC: 0xffffff,
  LV: 0xffffff,
  UUH: 0xffffff,
  UUS: 0xffffff,
  TS: 0xffffff,
  UUO: 0xffffff,
  OG: 0xffffff,
});
export type ElementSymbolColors = typeof ElementSymbolColors;

const DefaultElementSymbolColor = Color(0xffffff);
const Description = 'Assigns a color to every atom according to its chemical element.';

export const ElementSymbolColorThemeParams = {
  carbonColor: PD.MappedStatic(
    'chain-id',
    {
      'chain-id': PD.Group(ChainIdColorThemeParams),
      'unit-index': PD.Group(UnitIndexColorThemeParams, { label: 'Chain Instance' }),
      'entity-id': PD.Group(EntityIdColorThemeParams),
      'entity-source': PD.Group(EntitySourceColorThemeParams),
      'operator-name': PD.Group(OperatorNameColorThemeParams),
      'model-index': PD.Group(ModelIndexColorThemeParams),
      'structure-index': PD.Group(StructureIndexColorThemeParams),
      'trajectory-index': PD.Group(TrajectoryIndexColorThemeParams),
      uniform: PD.Group(UniformColorThemeParams),
      'element-symbol': PD.EmptyGroup(),
    },
    { description: 'Use chain-id coloring for carbon atoms.' },
  ),
  saturation: PD.Numeric(0, { min: -6, max: 6, step: 0.1 }),
  lightness: PD.Numeric(0.2, { min: -6, max: 6, step: 0.1 }),
  colors: PD.MappedStatic('default', {
    default: PD.EmptyGroup(),
    custom: PD.Group(getColorMapParams(ElementSymbolColors)),
  }),
};
export type ElementSymbolColorThemeParams = typeof ElementSymbolColorThemeParams;
export function getElementSymbolColorThemeParams(ctx: ThemeDataContext) {
  return PD.clone(ElementSymbolColorThemeParams);
}

type ElementSymbolColorThemeProps = PD.Values<ElementSymbolColorThemeParams>;

export function elementSymbolColor(colorMap: ElementSymbolColors, element: ElementSymbol): Color {
  const c = colorMap[element as keyof ElementSymbolColors];
  return c === undefined ? DefaultElementSymbolColor : c;
}

function getCarbonTheme(ctx: ThemeDataContext, props: ElementSymbolColorThemeProps['carbonColor']) {
  switch (props.name) {
    case 'chain-id':
      return ChainIdColorTheme(ctx, props.params);
    case 'unit-index':
      return UnitIndexColorTheme(ctx, props.params);
    case 'entity-id':
      return EntityIdColorTheme(ctx, props.params);
    case 'entity-source':
      return EntitySourceColorTheme(ctx, props.params);
    case 'operator-name':
      return OperatorNameColorTheme(ctx, props.params);
    case 'model-index':
      return ModelIndexColorTheme(ctx, props.params);
    case 'structure-index':
      return StructureIndexColorTheme(ctx, props.params);
    case 'trajectory-index':
      return TrajectoryIndexColorTheme(ctx, props.params);
    case 'uniform':
      return UniformColorTheme(ctx, props.params);
    case 'element-symbol':
      return undefined;
    default:
      assertUnreachable(props);
  }
}

export function ElementSymbolColorTheme(
  ctx: ThemeDataContext,
  props: PD.Values<ElementSymbolColorThemeParams>,
): ColorTheme<ElementSymbolColorThemeParams> {
  const colorMap = getAdjustedColorMap(
    props.colors.name === 'default' ? ElementSymbolColors : props.colors.params,
    props.saturation,
    props.lightness,
  );

  const carbonTheme = getCarbonTheme(ctx, props.carbonColor);
  const carbonColor = carbonTheme?.color;
  const contextHash = carbonTheme?.contextHash ?? -1;

  function elementColor(element: ElementSymbol, location: Location) {
    return carbonColor && element === 'C' ? carbonColor(location, false) : elementSymbolColor(colorMap, element);
  }

  function color(location: Location): Color {
    if (StructureElement.Location.is(location)) {
      if (Unit.isAtomic(location.unit)) {
        const { type_symbol } = location.unit.model.atomicHierarchy.atoms;
        return elementColor(type_symbol.value(location.element), location);
      }
    } else if (Bond.isLocation(location)) {
      if (Unit.isAtomic(location.aUnit)) {
        const { type_symbol } = location.aUnit.model.atomicHierarchy.atoms;
        const element = type_symbol.value(location.aUnit.elements[location.aIndex]);
        return elementColor(element, location);
      }
    }
    return DefaultElementSymbolColor;
  }

  const granularity =
    props.carbonColor.name === 'operator-name' || props.carbonColor.name === 'unit-index' ? 'groupInstance' : 'group';

  return {
    factory: ElementSymbolColorTheme,
    granularity,
    preferSmoothing: true,
    color,
    props,
    contextHash,
    description: Description,
    legend: TableLegend(
      Object.keys(colorMap).map((name) => {
        return [name, (colorMap as any)[name] as Color] as [string, Color];
      }),
    ),
  };
}

export const ElementSymbolColorThemeProvider: ColorTheme.Provider<ElementSymbolColorThemeParams, 'element-symbol'> = {
  name: 'element-symbol',
  label: 'Element Symbol',
  category: ColorThemeCategory.Atom,
  factory: ElementSymbolColorTheme,
  getParams: getElementSymbolColorThemeParams,
  defaultValues: PD.getDefaultValues(ElementSymbolColorThemeParams),
  isApplicable: (ctx: ThemeDataContext) => !!ctx.structure,
};
