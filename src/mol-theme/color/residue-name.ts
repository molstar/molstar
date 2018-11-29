/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Color, ColorMap } from 'mol-util/color';
import { StructureElement, Unit, Link, ElementIndex } from 'mol-model/structure';
import { Location } from 'mol-model/location';
import { ColorTheme } from '../color';
import { ParamDefinition as PD } from 'mol-util/param-definition'
import { ThemeDataContext } from '../theme';
import { TableLegend } from 'mol-util/color/tables';

// protein colors from Jmol http://jmol.sourceforge.net/jscolors/
const ResidueNameColors = ColorMap({
    // standard amino acids
    'ALA': 0x8CFF8C,
    'ARG': 0x00007C,
    'ASN': 0xFF7C70,
    'ASP': 0xA00042,
    'CYS': 0xFFFF70,
    'GLN': 0xFF4C4C,
    'GLU': 0x660000,
    'GLY': 0xEEEEEE,
    'HIS': 0x7070FF,
    'ILE': 0x004C00,
    'LEU': 0x455E45,
    'LYS': 0x4747B8,
    'MET': 0xB8A042,
    'PHE': 0x534C52,
    'PRO': 0x525252,
    'SER': 0xFF7042,
    'THR': 0xB84C00,
    'TRP': 0x4F4600,
    'TYR': 0x8C704C,
    'VAL': 0xFF8CFF,

    // rna bases
    'A': 0xDC143C,  // Crimson Red
    'G': 0x32CD32,  // Lime Green
    'I': 0x9ACD32,  // Yellow Green
    'C': 0xFFD700,  // Gold Yellow
    'T': 0x4169E1,  // Royal Blue
    'U': 0x40E0D0,  // Turquoise Cyan

    // dna bases
    'DA': 0xDC143C,
    'DG': 0x32CD32,
    'DI': 0x9ACD32,
    'DC': 0xFFD700,
    'DT': 0x4169E1,
    'DU': 0x40E0D0,

    // peptide bases
    'APN': 0xDC143C,
    'GPN': 0x32CD32,
    'CPN': 0xFFD700,
    'TPN': 0x4169E1,
})

const DefaultResidueNameColor = Color(0xFF00FF)
const Description = 'Assigns a color to every residue according to its name.'

export const ResidueNameColorThemeParams = {}
export type ResidueNameColorThemeParams = typeof ResidueNameColorThemeParams
export function getResidueNameColorThemeParams(ctx: ThemeDataContext) {
    return ResidueNameColorThemeParams // TODO return copy
}

export function residueNameColor(residueName: string): Color {
    const c = (ResidueNameColors as { [k: string]: Color })[residueName];
    return c === undefined ? DefaultResidueNameColor : c
}

function getAtomicCompId(unit: Unit.Atomic, element: ElementIndex) {
    const { modifiedResidues } = unit.model.properties
    const compId = unit.model.atomicHierarchy.residues.auth_comp_id.value(unit.residueIndex[element])
    const parentId = modifiedResidues.parentId.get(compId)
    return parentId === undefined ? compId : parentId
}

function getCoarseCompId(unit: Unit.Spheres | Unit.Gaussians, element: ElementIndex) {
    const seqIdBegin = unit.coarseElements.seq_id_begin.value(element)
    const seqIdEnd = unit.coarseElements.seq_id_end.value(element)
    if (seqIdBegin === seqIdEnd) {
        const { modifiedResidues } = unit.model.properties
        const entityKey = unit.coarseElements.entityKey[element]
        const seq = unit.model.sequence.byEntityKey[entityKey]
        let compId = seq.compId.value(seqIdBegin - 1) // 1-indexed
        const parentId = modifiedResidues.parentId.get(compId)
        return parentId === undefined ? compId : parentId
    }
}

export function ResidueNameColorTheme(ctx: ThemeDataContext, props: PD.Values<ResidueNameColorThemeParams>): ColorTheme<ResidueNameColorThemeParams> {
    function color(location: Location): Color {
        if (StructureElement.isLocation(location)) {
            if (Unit.isAtomic(location.unit)) {
                return residueNameColor(getAtomicCompId(location.unit, location.element))
            } else {
                const compId = getCoarseCompId(location.unit, location.element)
                if (compId) return residueNameColor(compId)
            }
        } else if (Link.isLocation(location)) {
            if (Unit.isAtomic(location.aUnit)) {
                return residueNameColor(getAtomicCompId(location.aUnit, location.aUnit.elements[location.aIndex]))
            } else {
                const compId = getCoarseCompId(location.aUnit, location.aUnit.elements[location.aIndex])
                if (compId) return residueNameColor(compId)
            }
        }
        return DefaultResidueNameColor
    }

    return {
        factory: ResidueNameColorTheme,
        granularity: 'group',
        color,
        props,
        description: Description,
        legend: TableLegend(Object.keys(ResidueNameColors).map(name => {
            return [name, (ResidueNameColors as any)[name] as Color] as [string, Color]
        }).concat([[ 'Unknown', DefaultResidueNameColor ]]))
    }
}

export const ResidueNameColorThemeProvider: ColorTheme.Provider<ResidueNameColorThemeParams> = {
    label: 'Residue Name',
    factory: ResidueNameColorTheme,
    getParams: getResidueNameColorThemeParams,
    defaultValues: PD.getDefaultValues(ResidueNameColorThemeParams),
    isApplicable: (ctx: ThemeDataContext) => !!ctx.structure
}