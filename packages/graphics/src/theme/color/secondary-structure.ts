/**
 * Copyright (c) 2018-2022 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Color, ColorMap } from '@molstar/core/util/color';
import { StructureElement, Unit, Bond, type ElementIndex } from '@molstar/model/model/structure';
import type { Location } from '@molstar/model/model/location';
import type { ColorTheme } from '../color.js';
import { SecondaryStructureType, MoleculeType } from '@molstar/model/model/structure/model/types';
import { getElementMoleculeType } from '@molstar/model/model/structure/util';
import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import type { ThemeDataContext } from '../theme.js';
import { TableLegend } from '@molstar/core/util/legend';
import {
  SecondaryStructureProvider,
  type SecondaryStructureValue,
} from '@molstar/model/props/computed/secondary-structure';
import { getAdjustedColorMap } from '@molstar/core/util/color/color';
import { getColorMapParams } from '@molstar/core/util/color/params';
import type { CustomProperty } from '@molstar/model/props/common/custom-property';
import { hash2 } from '@molstar/core/data/util/hash-functions';
import { ColorThemeCategory } from './categories.js';

// from Jmol http://jmol.sourceforge.net/jscolors/ (shapely)
export const SecondaryStructureColors = ColorMap({
  alphaHelix: 0xff0080,
  threeTenHelix: 0xa00080,
  piHelix: 0x600080,
  betaTurn: 0x6080ff,
  betaStrand: 0xffc800,
  coil: 0xffffff,
  bend: 0x66d8c9 /* biting original color used 0x00FF00 */,
  turn: 0x00b266,

  dna: 0xae00fe,
  rna: 0xfd0162,

  carbohydrate: 0xa6a6fa,
});
export type SecondaryStructureColors = typeof SecondaryStructureColors;

const DefaultSecondaryStructureColor = Color(0x808080);
const Description = 'Assigns a color based on the type of secondary structure and basic molecule type.';

export const SecondaryStructureColorThemeParams = {
  saturation: PD.Numeric(-1, { min: -6, max: 6, step: 0.1 }),
  lightness: PD.Numeric(0, { min: -6, max: 6, step: 0.1 }),
  colors: PD.MappedStatic('default', {
    default: PD.EmptyGroup(),
    custom: PD.Group(getColorMapParams(SecondaryStructureColors)),
  }),
};
export type SecondaryStructureColorThemeParams = typeof SecondaryStructureColorThemeParams;
export function getSecondaryStructureColorThemeParams(ctx: ThemeDataContext) {
  return SecondaryStructureColorThemeParams; // TODO return copy
}

export function secondaryStructureColor(
  colorMap: SecondaryStructureColors,
  unit: Unit,
  element: ElementIndex,
  computedSecondaryStructure?: SecondaryStructureValue,
): Color {
  let secStrucType = SecondaryStructureType.create(SecondaryStructureType.Flag.None);
  if (computedSecondaryStructure && Unit.isAtomic(unit)) {
    const secondaryStructure = computedSecondaryStructure.get(unit.invariantId);
    if (secondaryStructure)
      secStrucType = secondaryStructure.type[secondaryStructure.getIndex(unit.residueIndex[element])];
  }

  if (SecondaryStructureType.is(secStrucType, SecondaryStructureType.Flag.Helix)) {
    if (SecondaryStructureType.is(secStrucType, SecondaryStructureType.Flag.Helix3Ten)) {
      return colorMap.threeTenHelix;
    } else if (SecondaryStructureType.is(secStrucType, SecondaryStructureType.Flag.HelixPi)) {
      return colorMap.piHelix;
    }
    return colorMap.alphaHelix;
  } else if (SecondaryStructureType.is(secStrucType, SecondaryStructureType.Flag.Beta)) {
    return colorMap.betaStrand;
  } else if (SecondaryStructureType.is(secStrucType, SecondaryStructureType.Flag.Bend)) {
    return colorMap.bend;
  } else if (SecondaryStructureType.is(secStrucType, SecondaryStructureType.Flag.Turn)) {
    return colorMap.turn;
  } else {
    const moleculeType = getElementMoleculeType(unit, element);
    if (moleculeType === MoleculeType.DNA) {
      return colorMap.dna;
    } else if (moleculeType === MoleculeType.RNA) {
      return colorMap.rna;
    } else if (moleculeType === MoleculeType.Saccharide) {
      return colorMap.carbohydrate;
    } else if (moleculeType === MoleculeType.Protein) {
      return colorMap.coil;
    }
  }
  return DefaultSecondaryStructureColor;
}

export function SecondaryStructureColorTheme(
  ctx: ThemeDataContext,
  props: PD.Values<SecondaryStructureColorThemeParams>,
): ColorTheme<SecondaryStructureColorThemeParams> {
  const computedSecondaryStructure = ctx.structure && SecondaryStructureProvider.get(ctx.structure);
  const contextHash = computedSecondaryStructure
    ? hash2(computedSecondaryStructure.id, computedSecondaryStructure.version)
    : -1;

  const colorMap = getAdjustedColorMap(
    props.colors.name === 'default' ? SecondaryStructureColors : props.colors.params,
    props.saturation,
    props.lightness,
  );

  function color(location: Location): Color {
    if (StructureElement.Location.is(location)) {
      return secondaryStructureColor(colorMap, location.unit, location.element, computedSecondaryStructure?.value);
    } else if (Bond.isLocation(location)) {
      return secondaryStructureColor(
        colorMap,
        location.aUnit,
        location.aUnit.elements[location.aIndex],
        computedSecondaryStructure?.value,
      );
    }
    return DefaultSecondaryStructureColor;
  }

  return {
    factory: SecondaryStructureColorTheme,
    granularity: 'group',
    preferSmoothing: true,
    color,
    props,
    contextHash,
    description: Description,
    legend: TableLegend(
      Object.keys(colorMap)
        .map((name) => {
          return [name, (colorMap as any)[name] as Color] as [string, Color];
        })
        .concat([['Other', DefaultSecondaryStructureColor]]),
    ),
  };
}

export const SecondaryStructureColorThemeProvider: ColorTheme.Provider<
  SecondaryStructureColorThemeParams,
  'secondary-structure'
> = {
  name: 'secondary-structure',
  label: 'Secondary Structure',
  category: ColorThemeCategory.Residue,
  factory: SecondaryStructureColorTheme,
  getParams: getSecondaryStructureColorThemeParams,
  defaultValues: PD.getDefaultValues(SecondaryStructureColorThemeParams),
  isApplicable: (ctx: ThemeDataContext) => !!ctx.structure,
  ensureCustomProperties: {
    attach: (ctx: CustomProperty.Context, data: ThemeDataContext) =>
      data.structure ? SecondaryStructureProvider.attach(ctx, data.structure, void 0, true) : Promise.resolve(),
    detach: (data) => data.structure && SecondaryStructureProvider.ref(data.structure, false),
  },
};
