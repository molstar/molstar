/**
 * Copyright (c) 2018-2025 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author Lukáš Polák <admin@lukaspolak.cz>
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

// Colors for charged residues (by-name)
export const ChargedResidueColors = ColorMap({
  // standard amino acids (charged)
  ARG: 0x0000ff,
  ASP: 0xff0000,
  GLU: 0xff0000,
  HIS: 0x33c3f9,
  LYS: 0x0000ff,

  // standard amino acids (uncharged)
  ALA: 0xffffff,
  ASN: 0xffffff,
  CYS: 0xffffff,
  GLN: 0xffffff,
  GLY: 0xffffff,
  ILE: 0xffffff,
  LEU: 0xffffff,
  MET: 0xffffff,
  PHE: 0xffffff,
  PRO: 0xffffff,
  SER: 0xffffff,
  THR: 0xffffff,
  TRP: 0xffffff,
  TYR: 0xffffff,
  VAL: 0xffffff,

  // common from CCD
  MSE: 0xffffff,
  SEP: 0xffffff,
  TPO: 0xffffff,
  PTR: 0xffffff,
  PCA: 0xffffff,
  HYP: 0xffffff,

  // charmm ff
  HSD: 0xffffff,
  HSE: 0xffffff,
  HSP: 0x0000ff,
  LSN: 0xffffff,
  ASPP: 0xffffff,
  GLUP: 0xffffff,

  // amber ff
  HID: 0xffffff,
  HIE: 0xffffff,
  HIP: 0x0000ff,
  LYN: 0xffffff,
  ASH: 0xffffff,
  GLH: 0xffffff,

  // rna bases
  A: 0xffffff,
  G: 0xffffff,
  I: 0xffffff,
  C: 0xffffff,
  T: 0xffffff,
  U: 0xffffff,

  // dna bases
  DA: 0xffffff,
  DG: 0xffffff,
  DI: 0xffffff,
  DC: 0xffffff,
  DT: 0xffffff,
  DU: 0xffffff,

  // peptide bases
  APN: 0xffffff,
  GPN: 0xffffff,
  CPN: 0xffffff,
  TPN: 0xffffff,
});
export type ChargedResidueColors = typeof ChargedResidueColors;

const DefaultResidueChargeColor = Color(0xff00ff);
const Description = 'Assigns a color to every residue based on its charge state.';

export const ResidueChargeColorThemeParams = {
  method: PD.MappedStatic('by-name', {
    'by-name': PD.Group(
      {
        saturation: PD.Numeric(0, { min: -6, max: 6, step: 0.1 }),
        lightness: PD.Numeric(0, { min: -6, max: 6, step: 0.1 }),
        colors: PD.MappedStatic('default', {
          default: PD.EmptyGroup(),
          custom: PD.Group(getColorMapParams(ChargedResidueColors)),
        }),
      },
      { isFlat: true },
    ),
  }),
};
export type ResidueChargeColorThemeParams = typeof ResidueChargeColorThemeParams;
export function getResidueChargeColorThemeParams(ctx: ThemeDataContext) {
  return PD.clone(ResidueChargeColorThemeParams);
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

export function residueChargeColor(colorMap: ColorMap<Record<string, Color>>, residueName: string): Color {
  const c = colorMap[residueName];
  return c === undefined ? DefaultResidueChargeColor : c;
}

export function ResidueChargeColorTheme(
  ctx: ThemeDataContext,
  props: PD.Values<ResidueChargeColorThemeParams>,
): ColorTheme<ResidueChargeColorThemeParams> {
  const { saturation, lightness, colors } = props.method.params;
  const colorMap = getAdjustedColorMap(
    props.method.params.colors.name === 'default' ? ChargedResidueColors : colors.params,
    saturation,
    lightness,
  );

  function color(location: Location): Color {
    if (StructureElement.Location.is(location)) {
      if (Unit.isAtomic(location.unit)) {
        const compId = getAtomicCompId(location.unit, location.element);
        return residueChargeColor(colorMap, compId);
      } else {
        const compId = getCoarseCompId(location.unit, location.element);
        if (compId) return residueChargeColor(colorMap, compId);
      }
    } else if (Bond.isLocation(location)) {
      if (Unit.isAtomic(location.aUnit)) {
        const compId = getAtomicCompId(location.aUnit, location.aUnit.elements[location.aIndex]);
        return residueChargeColor(colorMap, compId);
      } else {
        const compId = getCoarseCompId(location.aUnit, location.aUnit.elements[location.aIndex]);
        if (compId) return residueChargeColor(colorMap, compId);
      }
    }
    return DefaultResidueChargeColor;
  }

  return {
    factory: ResidueChargeColorTheme,
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
        .concat([['Unknown', DefaultResidueChargeColor]]),
    ),
  };
}

export const ResidueChargeColorThemeProvider: ColorTheme.Provider<ResidueChargeColorThemeParams, 'residue-charge'> = {
  name: 'residue-charge',
  label: 'Residue Charge',
  category: ColorThemeCategory.Residue,
  factory: ResidueChargeColorTheme,
  getParams: getResidueChargeColorThemeParams,
  defaultValues: PD.getDefaultValues(ResidueChargeColorThemeParams),
  isApplicable: (ctx: ThemeDataContext) => !!ctx.structure,
};
