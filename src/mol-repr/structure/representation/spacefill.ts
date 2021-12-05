/**
 * Copyright (c) 2018-2021 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { ElementSphereVisual, ElementSphereParams } from '../visual/element-sphere';
import { UnitsRepresentation } from '../units-representation';
import { ParamDefinition as PD } from '../../../mol-util/param-definition';
import { StructureRepresentation, StructureRepresentationProvider, StructureRepresentationStateBuilder } from '../representation';
import { RepresentationParamsGetter, RepresentationContext, Representation } from '../../../mol-repr/representation';
import { ThemeRegistryContext } from '../../../mol-theme/theme';
import { Structure } from '../../../mol-model/structure';
import { BaseGeometry } from '../../../mol-geo/geometry/base';

const SpacefillVisuals = {
    'element-sphere': (ctx: RepresentationContext, getParams: RepresentationParamsGetter<Structure, ElementSphereParams>) => UnitsRepresentation('Sphere mesh', ctx, getParams, ElementSphereVisual),
};

export const SpacefillParams = {
    ...ElementSphereParams,
    bumpFrequency: PD.Numeric(1, { min: 0, max: 10, step: 0.1 }, BaseGeometry.ShadingCategory),
};
export type SpacefillParams = typeof SpacefillParams
export function getSpacefillParams(ctx: ThemeRegistryContext, structure: Structure) {
    const params = PD.clone(SpacefillParams);
    if (structure.isCoarseGrained) {
        params.sizeFactor.defaultValue = 2;
    }
    return params;
}

export type SpacefillRepresentation = StructureRepresentation<SpacefillParams>
export function SpacefillRepresentation(ctx: RepresentationContext, getParams: RepresentationParamsGetter<Structure, SpacefillParams>): SpacefillRepresentation {
    return Representation.createMulti('Spacefill', ctx, getParams, StructureRepresentationStateBuilder, SpacefillVisuals as unknown as Representation.Def<Structure, SpacefillParams>);
}

export const SpacefillRepresentationProvider = StructureRepresentationProvider({
    name: 'spacefill',
    label: 'Spacefill',
    description: 'Displays atomic/coarse elements as spheres.',
    factory: SpacefillRepresentation,
    getParams: getSpacefillParams,
    defaultValues: PD.getDefaultValues(SpacefillParams),
    defaultColorTheme: { name: 'element-symbol' },
    defaultSizeTheme: { name: 'physical' },
    isApplicable: (structure: Structure) => structure.elementCount > 0
});