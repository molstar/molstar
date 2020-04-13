/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { ElementPointVisual, ElementPointParams } from '../visual/element-point';
import { UnitsRepresentation } from '../units-representation';
import { ParamDefinition as PD } from '../../../mol-util/param-definition';
import { StructureRepresentation, StructureRepresentationProvider, StructureRepresentationStateBuilder } from '../representation';
import { Representation, RepresentationParamsGetter, RepresentationContext } from '../../../mol-repr/representation';
import { ThemeRegistryContext } from '../../../mol-theme/theme';
import { Structure } from '../../../mol-model/structure';

const PointVisuals = {
    'element-point': (ctx: RepresentationContext, getParams: RepresentationParamsGetter<Structure, ElementPointParams>) => UnitsRepresentation('Points', ctx, getParams, ElementPointVisual),
};

export const PointParams = {
    ...ElementPointParams,
};
export type PointParams = typeof PointParams
export function getPointParams(ctx: ThemeRegistryContext, structure: Structure) {
    return PD.clone(PointParams);
}

export type PointRepresentation = StructureRepresentation<PointParams>
export function PointRepresentation(ctx: RepresentationContext, getParams: RepresentationParamsGetter<Structure, PointParams>): PointRepresentation {
    return Representation.createMulti('Point', ctx, getParams, StructureRepresentationStateBuilder, PointVisuals as unknown as Representation.Def<Structure, PointParams>);
}

export const PointRepresentationProvider = StructureRepresentationProvider({
    name: 'point',
    label: 'Point',
    description: 'Displays elements (atoms, coarse spheres) as spheres.',
    factory: PointRepresentation,
    getParams: getPointParams,
    defaultValues: PD.getDefaultValues(PointParams),
    defaultColorTheme: { name: 'element-symbol' },
    defaultSizeTheme: { name: 'physical' },
    isApplicable: (structure: Structure) => structure.elementCount > 0
});