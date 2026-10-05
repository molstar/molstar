/**
 * Copyright (c) 2019-2024 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { PolymerTubeVisual, PolymerTubeParams } from '../visual/polymer-tube-mesh.js';
import { PolymerGapVisual, PolymerGapParams } from '../visual/polymer-gap-cylinder.js';
import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import { UnitsRepresentation } from '../units-representation.js';
import { type StructureRepresentation, StructureRepresentationProvider, StructureRepresentationStateBuilder } from '../representation.js';
import { Representation, type RepresentationParamsGetter, type RepresentationContext } from '../../representation.js';
import { type Structure, Unit } from '@molstar/model/model/structure';
import type { ThemeRegistryContext } from '@molstar/graphics/theme/theme';
import { BaseGeometry } from '@molstar/graphics/geo/geometry/base';

const PuttyVisuals = {
    'polymer-tube': (ctx: RepresentationContext, getParams: RepresentationParamsGetter<Structure, PolymerTubeParams>) => UnitsRepresentation('Polymer tube mesh', ctx, getParams, PolymerTubeVisual),
    'polymer-gap': (ctx: RepresentationContext, getParams: RepresentationParamsGetter<Structure, PolymerGapParams>) => UnitsRepresentation('Polymer gap cylinder', ctx, getParams, PolymerGapVisual),
};

export const PuttyParams = {
    ...PolymerTubeParams,
    ...PolymerGapParams,
    sizeFactor: PD.Numeric(0.2, { min: 0, max: 10, step: 0.01 }),
    visuals: PD.MultiSelect(['polymer-tube', 'polymer-gap'], PD.objectToOptions(PuttyVisuals)),
    bumpFrequency: PD.Numeric(2, { min: 0, max: 10, step: 0.1 }, BaseGeometry.ShadingCategory),
    density: PD.Numeric(0.1, { min: 0, max: 1, step: 0.01 }, BaseGeometry.ShadingCategory),
};
export type PuttyParams = typeof PuttyParams
export function getPuttyParams(ctx: ThemeRegistryContext, structure: Structure) {
    const params = PD.clone(PuttyParams);
    let hasNucleotides = false;
    let hasGaps = false;
    structure.units.forEach(u => {
        if (!hasNucleotides && Unit.isAtomic(u) && u.nucleotideElements.length) hasNucleotides = true;
        if (!hasGaps && u.gapElements.length) hasGaps = true;
    });
    params.visuals.defaultValue = ['polymer-tube'];
    if (hasGaps) params.visuals.defaultValue.push('polymer-gap');
    return params;
}

export type PuttyRepresentation = StructureRepresentation<PuttyParams>
export function PuttyRepresentation(ctx: RepresentationContext, getParams: RepresentationParamsGetter<Structure, PuttyParams>): PuttyRepresentation {
    return Representation.createMulti('Putty', ctx, getParams, StructureRepresentationStateBuilder, PuttyVisuals as unknown as Representation.Def<Structure, PuttyParams>);
}

export const PuttyRepresentationProvider = StructureRepresentationProvider({
    name: 'putty',
    label: 'Putty',
    description: 'Displays a tube smoothly following the trace atoms of polymers.',
    factory: PuttyRepresentation,
    getParams: getPuttyParams,
    defaultValues: PD.getDefaultValues(PuttyParams),
    defaultColorTheme: { name: 'chain-id' },
    defaultSizeTheme: { name: 'uncertainty' },
    isApplicable: (structure: Structure) => structure.polymerResidueCount > 0,
});