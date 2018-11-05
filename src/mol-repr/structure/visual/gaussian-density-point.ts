/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Unit, Structure } from 'mol-model/structure';
import { UnitsVisual } from '../index';
import { VisualUpdateState } from '../../util';
import { StructureElementIterator } from './util/element';
import { EmptyLoci } from 'mol-model/loci';
import { Vec3 } from 'mol-math/linear-algebra';
import { UnitsPointsVisual, UnitsPointsParams } from '../units-visual';
import { SizeThemeOptions, SizeThemeName } from 'mol-theme/size';
import { GaussianDensityProps, GaussianDensityParams } from 'mol-model/structure/structure/unit/gaussian-density';
import { paramDefaultValues, SelectParam, NumberParam, BooleanParam } from 'mol-util/parameter';
import { Points } from 'mol-geo/geometry/points/points';
import { PointsBuilder } from 'mol-geo/geometry/points/points-builder';
import { VisualContext } from 'mol-repr';

export const GaussianDensityPointParams = {
    ...UnitsPointsParams,
    ...GaussianDensityParams,
    sizeTheme: SelectParam<SizeThemeName>('Size Theme', '', 'uniform', SizeThemeOptions),
    sizeValue: NumberParam('Size Value', '', 1, 0, 20, 0.1),
    pointSizeAttenuation: BooleanParam('Point Size Attenuation', '', false),
}
export const DefaultGaussianDensityPointProps = paramDefaultValues(GaussianDensityPointParams)
export type GaussianDensityPointProps = typeof DefaultGaussianDensityPointProps

export async function createGaussianDensityPoint(ctx: VisualContext, unit: Unit, structure: Structure, props: GaussianDensityProps, points?: Points) {
    const { transform, field: { space, data } } = await unit.computeGaussianDensity(props, ctx.runtime, ctx.webgl)

    const { dimensions, get } = space
    const [ xn, yn, zn ] = dimensions

    const n = xn * yn * zn * 3
    const builder = PointsBuilder.create(n, n / 10, points)

    const p = Vec3.zero()
    let i = 0

    for (let x = 0; x < xn; ++x) {
        for (let y = 0; y < yn; ++y) {
            for (let z = 0; z < zn; ++z) {
                if (get(data, x, y, z) > 0.001) {
                    Vec3.set(p, x, y, z)
                    Vec3.transformMat4(p, p, transform)
                    builder.add(p[0], p[1], p[2], i)
                }
                if (i % 100000 === 0 && ctx.runtime.shouldUpdate) {
                    await ctx.runtime.update({ message: 'Creating density points', current: i, max: n });
                }
                ++i
            }
        }
    }
    return builder.getPoints()
}

export function GaussianDensityPointVisual(): UnitsVisual<GaussianDensityPointProps> {
    return UnitsPointsVisual<GaussianDensityPointProps>({
        defaultProps: DefaultGaussianDensityPointProps,
        createGeometry: createGaussianDensityPoint,
        createLocationIterator: StructureElementIterator.fromGroup,
        getLoci: () => EmptyLoci,
        mark: () => false,
        setUpdateState: (state: VisualUpdateState, newProps: GaussianDensityPointProps, currentProps: GaussianDensityPointProps) => {
            if (newProps.resolution !== currentProps.resolution) state.createGeometry = true
            if (newProps.radiusOffset !== currentProps.radiusOffset) state.createGeometry = true
            if (newProps.smoothness !== currentProps.smoothness) state.createGeometry = true
            if (newProps.useGpu !== currentProps.useGpu) state.createGeometry = true
            if (newProps.ignoreCache !== currentProps.ignoreCache) state.createGeometry = true
        }
    })
}