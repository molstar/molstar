/**
 * Copyright (c) 2018-2022 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Color, ColorMap } from '@molstar/core/util/color';
import { StructureElement, Unit, Bond, type ElementIndex } from '@molstar/model/model/structure';
import type { Location } from '@molstar/model/model/location';
import type { ColorTheme } from '../color.js';
import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import type { ThemeDataContext } from '../theme.js';
import { TableLegend } from '@molstar/core/util/legend';
import { getAdjustedColorMap } from '@molstar/core/util/color/color';
import { getColorMapParams } from '@molstar/core/util/color/params';
import { ColorThemeCategory } from './categories.js';

// protein colors from Jmol http://jmol.sourceforge.net/jscolors/
export const ResidueNameColors = ColorMap({
  // standard amino acids
  ALA: 0x8cff8c,
  ARG: 0x00007c,
  ASN: 0xff7c70,
  ASP: 0xa00042,
  CYS: 0xffff70,
  GLN: 0xff4c4c,
  GLU: 0x660000,
  GLY: 0xeeeeee,
  HIS: 0x7070ff,
  ILE: 0x004c00,
  LEU: 0x455e45,
  LYS: 0x4747b8,
  MET: 0xb8a042,
  PHE: 0x534c52,
  PRO: 0x525252,
  SER: 0xff7042,
  THR: 0xb84c00,
  TRP: 0x4f4600,
  TYR: 0x8c704c,
  VAL: 0xff8cff,

  // rna bases
  A: 0xdc143c, // Crimson Red
  G: 0x32cd32, // Lime Green
  I: 0x9acd32, // Yellow Green
  C: 0xffd700, // Gold Yellow
  T: 0x4169e1, // Royal Blue
  U: 0x40e0d0, // Turquoise Cyan

  // dna bases
  DA: 0xdc143c,
  DG: 0x32cd32,
  DI: 0x9acd32,
  DC: 0xffd700,
  DT: 0x4169e1,
  DU: 0x40e0d0,

  // peptide bases
  APN: 0xdc143c,
  GPN: 0x32cd32,
  CPN: 0xffd700,
  TPN: 0x4169e1,
});
export type ResidueNameColors = typeof ResidueNameColors;

const DefaultResidueNameColor = Color(0xff00ff);
const Description = 'Assigns a color to every residue according to its name.';

export const ResidueNameColorThemeParams = {
  saturation: PD.Numeric(0, { min: -6, max: 6, step: 0.1 }),
  lightness: PD.Numeric(1, { min: -6, max: 6, step: 0.1 }),
  colors: PD.MappedStatic('default', {
    default: PD.EmptyGroup(),
    custom: PD.Group(getColorMapParams(ResidueNameColors)),
  }),
};
export type ResidueNameColorThemeParams = typeof ResidueNameColorThemeParams;
export function getResidueNameColorThemeParams(ctx: ThemeDataContext) {
  return ResidueNameColorThemeParams; // TODO return copy
}

function getAtomicCompId(unit: Unit.Atomic, element: ElementIndex) {
  return unit.model.atomicHierarchy.atoms.label_comp_id.value(element);
}

function getCoarseCompId(unit: Unit.Spheres | Unit.Gaussians, element: ElementIndex) {
  const seqIdBegin = unit.coarseElements.seq_id_begin.value(element);
  const seqIdEnd = unit.coarseElements.seq_id_end.value(element);
  if (seqIdBegin === seqIdEnd) {
    const entityKey = unit.coarseElements.entityKey[element];
    const seq = unit.model.sequence.byEntityKey[entityKey].sequence;
    return seq.compId.value(seqIdBegin - 1); // 1-indexed
  }
}

export function residueNameColor(colorMap: ResidueNameColors, residueName: string): Color {
  const c = colorMap[residueName as keyof ResidueNameColors];
  return c === undefined ? DefaultResidueNameColor : c;
}

export function ResidueNameColorTheme(
  ctx: ThemeDataContext,
  props: PD.Values<ResidueNameColorThemeParams>,
): ColorTheme<ResidueNameColorThemeParams> {
  const colorMap = getAdjustedColorMap(
    props.colors.name === 'default' ? ResidueNameColors : props.colors.params,
    props.saturation,
    props.lightness,
  );

  function color(location: Location): Color {
    if (StructureElement.Location.is(location)) {
      if (Unit.isAtomic(location.unit)) {
        const compId = getAtomicCompId(location.unit, location.element);
        return residueNameColor(colorMap, compId);
      } else {
        const compId = getCoarseCompId(location.unit, location.element);
        if (compId) return residueNameColor(colorMap, compId);
      }
    } else if (Bond.isLocation(location)) {
      if (Unit.isAtomic(location.aUnit)) {
        const compId = getAtomicCompId(location.aUnit, location.aUnit.elements[location.aIndex]);
        return residueNameColor(colorMap, compId);
      } else {
        const compId = getCoarseCompId(location.aUnit, location.aUnit.elements[location.aIndex]);
        if (compId) return residueNameColor(colorMap, compId);
      }
    }
    return DefaultResidueNameColor;
  }

  return {
    factory: ResidueNameColorTheme,
    granularity: 'group',
    preferSmoothing: true,
    color,
    props,
    description: Description,
    legend: TableLegend(
      Object.keys(colorMap)
        .map((name) => {
          return [name, (colorMap as any)[name] as Color] as [string, Color];
        })
        .concat([['Unknown', DefaultResidueNameColor]]),
    ),
  };
}

export const ResidueNameColorThemeProvider: ColorTheme.Provider<ResidueNameColorThemeParams, 'residue-name'> = {
  name: 'residue-name',
  label: 'Residue Name',
  category: ColorThemeCategory.Residue,
  factory: ResidueNameColorTheme,
  getParams: getResidueNameColorThemeParams,
  defaultValues: PD.getDefaultValues(ResidueNameColorThemeParams),
  isApplicable: (ctx: ThemeDataContext) => !!ctx.structure,
};
