/**
 * Copyright (c) 2025 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { ParamDefinition as PD } from '../../../mol-util/param-definition';
import { ComplexRepresentation, StructureRepresentation, StructureRepresentationProvider, StructureRepresentationStateBuilder } from '../representation';
import { Representation, RepresentationParamsGetter, RepresentationContext } from '../../representation';
import { ThemeRegistryContext } from '../../../mol-theme/theme';
import { Structure } from '../../../mol-model/structure';
import { PlaneImageParams, PlaneImageVisual } from '../visual/plane-image';

const PlaneVisuals = {
    'plane-image': (ctx: RepresentationContext, getParams: RepresentationParamsGetter<Structure, PlaneImageParams>) => ComplexRepresentation('Plane image', ctx, getParams, PlaneImageVisual),
};

export const PlaneParams = {
    ...PlaneImageParams,
    visuals: PD.MultiSelect(['plane-image'], PD.objectToOptions(PlaneVisuals)),
};
export type PlaneParams = typeof PlaneParams
export function getPlaneParams(ctx: ThemeRegistryContext, structure: Structure) {
    return PlaneParams;
}

export type PlaneRepresentation = StructureRepresentation<PlaneParams>
export function PlaneRepresentation(ctx: RepresentationContext, getParams: RepresentationParamsGetter<Structure, PlaneParams>): PlaneRepresentation {
    return Representation.createMulti('Plane', ctx, getParams, StructureRepresentationStateBuilder, PlaneVisuals as unknown as Representation.Def<Structure, PlaneParams>);
}

export const PlaneRepresentationProvider = StructureRepresentationProvider({
    name: 'plane',
    label: 'Plane',
    description: 'Displays planes.',
    factory: PlaneRepresentation,
    getParams: getPlaneParams,
    defaultValues: PD.getDefaultValues(PlaneParams),
    defaultColorTheme: { name: 'element-symbol' },
    defaultSizeTheme: { name: 'physical' },
    isApplicable: (structure: Structure) => structure.elementCount > 0,
    getData: (structure: Structure, props: PD.Values<PlaneParams>) => {
        return props.includeParent ? structure.asParent() : structure;
    },
    mustRecreate: (oldProps: PD.Values<PlaneParams>, newProps: PD.Values<PlaneParams>) => {
        return oldProps.includeParent !== newProps.includeParent;
    }
});