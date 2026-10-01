/**
 * Copyright (c) 2018-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { ElementSphereVisual, ElementSphereParams, StructureElementSphereVisual } from '../visual/element-sphere.js';
import { UnitsRepresentation } from '../units-representation.js';
import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import { ComplexRepresentation, type StructureRepresentation, StructureRepresentationProvider, StructureRepresentationStateBuilder } from '../representation.js';
import { type RepresentationParamsGetter, type RepresentationContext, Representation } from '../../representation.js';
import type { ThemeRegistryContext } from '@molstar/graphics/theme/theme';
import type { Structure } from '@molstar/model/model/structure';
import { BaseGeometry } from '@molstar/graphics/geo/geometry/base';

const SpacefillVisuals = {
    'element-sphere': (ctx: RepresentationContext, getParams: RepresentationParamsGetter<Structure, ElementSphereParams>) => UnitsRepresentation('Sphere mesh/impostor', ctx, getParams, ElementSphereVisual),
    'structure-element-sphere': (ctx: RepresentationContext, getParams: RepresentationParamsGetter<Structure, ElementSphereParams>) => ComplexRepresentation('Structure sphere mesh/impostor', ctx, getParams, StructureElementSphereVisual),
};

export const SpacefillParams = {
    ...ElementSphereParams,
    bumpFrequency: PD.Numeric(1, { min: 0, max: 10, step: 0.1 }, BaseGeometry.ShadingCategory),
    density: PD.Numeric(0.5, { min: 0, max: 1, step: 0.01 }, BaseGeometry.ShadingCategory),
    visuals: PD.MultiSelect(['element-sphere'], PD.objectToOptions(SpacefillVisuals)),
};
export type SpacefillParams = typeof SpacefillParams

let CoarseGrainedSpacefillParams: SpacefillParams;
export function getSpacefillParams(ctx: ThemeRegistryContext, structure: Structure) {
    let params = SpacefillParams;
    if (structure.isCoarseGrained) {
        if (!CoarseGrainedSpacefillParams) {
            CoarseGrainedSpacefillParams = PD.clone(SpacefillParams);
            CoarseGrainedSpacefillParams.sizeFactor.defaultValue = 2;
        }
        params = CoarseGrainedSpacefillParams;
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