/**
 * Copyright (c) 2018-2024 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { PolymerBackboneCylinderVisual, PolymerBackboneCylinderParams } from '../visual/polymer-backbone-cylinder';
import { ParamDefinition as PD } from '../../../mol-util/param-definition';
import { UnitsRepresentation } from '../units-representation';
import { StructureRepresentation, StructureRepresentationProvider, StructureRepresentationStateBuilder } from '../representation';
import { Representation, RepresentationContext, RepresentationParamsGetter } from '../../../mol-repr/representation';
import { ThemeRegistryContext } from '../../../mol-theme/theme';
import { Structure } from '../../../mol-model/structure';
import { PolymerBackboneSphereParams, PolymerBackboneSphereVisual } from '../visual/polymer-backbone-sphere';
import { PolymerGapParams, PolymerGapVisual } from '../visual/polymer-gap-cylinder';
import { BaseGeometry } from '../../../mol-geo/geometry/base';

const BackboneVisuals = {
    'polymer-backbone-cylinder': (ctx: RepresentationContext, getParams: RepresentationParamsGetter<Structure, PolymerBackboneCylinderParams>) => UnitsRepresentation('Polymer backbone cylinder', ctx, getParams, PolymerBackboneCylinderVisual),
    'polymer-backbone-sphere': (ctx: RepresentationContext, getParams: RepresentationParamsGetter<Structure, PolymerBackboneSphereParams>) => UnitsRepresentation('Polymer backbone sphere', ctx, getParams, PolymerBackboneSphereVisual),
    'polymer-gap': (ctx: RepresentationContext, getParams: RepresentationParamsGetter<Structure, PolymerGapParams>) => UnitsRepresentation('Polymer gap cylinder', ctx, getParams, PolymerGapVisual),
};

export const BackboneParams = {
    ...PolymerBackboneSphereParams,
    ...PolymerBackboneCylinderParams,
    ...PolymerGapParams,
    sizeAspectRatio: PD.Numeric(1, { min: 0.1, max: 3, step: 0.1 }),
    visuals: PD.MultiSelect(['polymer-backbone-cylinder', 'polymer-backbone-sphere', 'polymer-gap'], PD.objectToOptions(BackboneVisuals)),
    bumpFrequency: PD.Numeric(0, { min: 0, max: 10, step: 0.1 }, BaseGeometry.ShadingCategory),
    density: PD.Numeric(0.1, { min: 0, max: 1, step: 0.01 }, BaseGeometry.ShadingCategory),
    colorMode: PD.Select('default', PD.arrayToOptions(['default', 'interpolate'] as const), { ...BaseGeometry.ShadingCategory, isHidden: true }),
};
export type BackboneParams = typeof BackboneParams
export function getBackboneParams(ctx: ThemeRegistryContext, structure: Structure) {
    const params = PD.clone(BackboneParams);
    let hasGaps = false;
    structure.units.forEach(u => {
        if (!hasGaps && u.gapElements.length) hasGaps = true;
    });
    params.visuals.defaultValue = ['polymer-backbone-cylinder', 'polymer-backbone-sphere'];
    if (hasGaps) params.visuals.defaultValue.push('polymer-gap');
    return params;
}

export type BackboneRepresentation = StructureRepresentation<BackboneParams>
export function BackboneRepresentation(ctx: RepresentationContext, getParams: RepresentationParamsGetter<Structure, BackboneParams>): BackboneRepresentation {
    return Representation.createMulti('Backbone', ctx, getParams, StructureRepresentationStateBuilder, BackboneVisuals as unknown as Representation.Def<Structure, BackboneParams>);
}

export const BackboneRepresentationProvider = StructureRepresentationProvider({
    name: 'backbone',
    label: 'Backbone',
    description: 'Displays polymer backbone with cylinders and spheres.',
    factory: BackboneRepresentation,
    getParams: getBackboneParams,
    defaultValues: PD.getDefaultValues(BackboneParams),
    defaultColorTheme: { name: 'chain-id' },
    defaultSizeTheme: { name: 'uniform' },
    isApplicable: (structure: Structure) => structure.polymerResidueCount > 0,
});