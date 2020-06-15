/**
 * Copyright (c) 2019-2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { PolymerTubeVisual,  PolymerTubeParams } from '../visual/polymer-tube-mesh';
import { PolymerGapVisual, PolymerGapParams } from '../visual/polymer-gap-cylinder';
import { ParamDefinition as PD } from '../../../mol-util/param-definition';
import { UnitsRepresentation } from '../units-representation';
import { StructureRepresentation, StructureRepresentationProvider, StructureRepresentationStateBuilder } from '../representation';
import { Representation, RepresentationParamsGetter, RepresentationContext } from '../../../mol-repr/representation';
import { Structure, Unit } from '../../../mol-model/structure';
import { ThemeRegistryContext } from '../../../mol-theme/theme';
import { CustomProperty } from '../../../mol-model-props/common/custom-property';
import { SecondaryStructureProvider } from '../../../mol-model-props/computed/secondary-structure';

const PuttyVisuals = {
    'polymer-tube': (ctx: RepresentationContext, getParams: RepresentationParamsGetter<Structure, PolymerTubeParams>) => UnitsRepresentation('Polymer tube mesh', ctx, getParams, PolymerTubeVisual),
    'polymer-gap': (ctx: RepresentationContext, getParams: RepresentationParamsGetter<Structure, PolymerGapParams>) => UnitsRepresentation('Polymer gap cylinder', ctx, getParams, PolymerGapVisual),
};

export const PuttyParams = {
    ...PolymerTubeParams,
    ...PolymerGapParams,
    sizeFactor: PD.Numeric(0.2, { min: 0, max: 10, step: 0.01 }),
    visuals: PD.MultiSelect(['polymer-tube', 'polymer-gap'], PD.objectToOptions(PuttyVisuals)),
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
    ensureCustomProperties: {
        attach: (ctx: CustomProperty.Context, structure: Structure) => SecondaryStructureProvider.attach(ctx, structure, void 0, true),
        detach: (data) => SecondaryStructureProvider.ref(data, false)
    }
});