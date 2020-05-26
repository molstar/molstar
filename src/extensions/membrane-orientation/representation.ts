/**
 * Copyright (c) 2018-2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author Sebastian Bittrich <sebastian.bittrich@rcsb.org>
 */

import { ParamDefinition as PD } from '../../mol-util/param-definition';
import { Vec3 } from '../../mol-math/linear-algebra';
import { Representation, RepresentationContext, RepresentationParamsGetter } from '../../mol-repr/representation';
import { Structure } from '../../mol-model/structure';
import { Spheres } from '../../mol-geo/geometry/spheres/spheres';
import { SpheresBuilder } from '../../mol-geo/geometry/spheres/spheres-builder';
import { MembraneOrientationProvider } from '../../mol-model-props/computed/membrane-orientation';
import { StructureRepresentationProvider, StructureRepresentation, StructureRepresentationStateBuilder } from '../../mol-repr/structure/representation';
import { MembraneOrientation } from '../../mol-model/structure/model/properties/membrane-orientation';
import { ThemeRegistryContext } from '../../mol-theme/theme';
import { ShapeRepresentation } from '../../mol-repr/shape/representation';
import { Shape } from '../../mol-model/shape';
import { ColorNames } from '../../mol-util/color/names';
import { RuntimeContext } from '../../mol-task';

const MembraneOrientationParams = {
    ...Spheres.Params,
    color: PD.Color(ColorNames.lightgrey),
    size: PD.Numeric(1, { min: 0.1, max: 10, step: 0.1 }, { description: 'Size of spheres that represent membrane planes' }),
    density: PD.Numeric(1, { min: 0.25, max: 10, step: 0.25 }, { description: 'Distance between spheres'})
};
export type MembraneOrientationParams = typeof MembraneOrientationParams
export type MembraneOrientationProps = PD.Values<MembraneOrientationParams>

export function getMembraneOrientationParams(ctx: ThemeRegistryContext, structure: Structure) {
    return PD.clone(MembraneOrientationParams);
}

export type MembraneOrientationRepresentation = StructureRepresentation<MembraneOrientationParams>
export function MembraneOrientationRepresentation(ctx: RepresentationContext, getParams: RepresentationParamsGetter<Structure, MembraneOrientationParams>): MembraneOrientationRepresentation {
    return Representation.createMulti('Membrane Orientation', ctx, getParams, StructureRepresentationStateBuilder, MembraneOrientationVisuals as unknown as Representation.Def<Structure, MembraneOrientationParams>);
}

const MembraneOrientationVisuals = {
    'membrane-orientation': (ctx: RepresentationContext, getParams: RepresentationParamsGetter<MembraneOrientation, MembraneOrientationParams>) => ShapeRepresentation(getMembraneSpheres, Spheres.Utils),
};


export const MembraneOrientationRepresentationProvider = StructureRepresentationProvider({
    name: 'membrane-orientation',
    label: 'Membrane Orientation',
    description: 'Displays a grid of points representing membrane layers.',
    factory: MembraneOrientationRepresentation,
    getParams: getMembraneOrientationParams,
    defaultValues: PD.getDefaultValues(MembraneOrientationParams),
    defaultColorTheme: { name: 'uniform' },
    defaultSizeTheme: { name: 'uniform' },
    isApplicable: (structure: Structure) => structure.elementCount > 0
});

function getMembraneSpheres(ctx: RuntimeContext, data: Structure, props: MembraneOrientationProps) {
    const { density } = props;
    const { radius, p1, p2, normal } = MembraneOrientationProvider.get(data).value!;;

    const spheresBuilder = SpheresBuilder.create();
    createMembraneLayer(spheresBuilder, p1, normal, density, radius);
    createMembraneLayer(spheresBuilder, p2, normal, density, radius);
    return Shape.create(name, data, spheresBuilder.getSpheres(), () => props.color, () => props.size, () => 'Membrane Boundaries');
}

function createMembraneLayer(spheresBuilder: SpheresBuilder, point: Vec3, normalVector: Vec3, density: number, radius: number) {
    const d = -Vec3.dot(normalVector, point);
    const rep = Vec3();
    for (let i = -1000, il = 1000; i < il; i += density) {
        for (let j = -1000, jl = 1000; j < jl; j += density) {
            Vec3.set(rep, i, j, -(d + i * normalVector[0] + j * normalVector[1]) / normalVector[2]);
            if (Vec3.squaredDistance(rep, point) < radius) {
                spheresBuilder.add(rep[0], rep[1], rep[2], 0);
            }
        }
    }
}