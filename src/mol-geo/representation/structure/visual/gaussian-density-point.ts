/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Unit, Structure } from 'mol-model/structure';
import { RuntimeContext } from 'mol-task'
import { UnitsVisual, VisualUpdateState } from '..';
import { StructureElementIterator } from './util/element';
import { EmptyLoci } from 'mol-model/loci';
import { Vec3 } from 'mol-math/linear-algebra';
import { SizeThemeProps } from 'mol-view/theme/size';
import { UnitsPointVisual, DefaultUnitsPointProps } from '../units-visual';
import { computeGaussianDensity } from './util/gaussian';
import { Color } from 'mol-util/color';
import { ColorThemeProps } from 'mol-view/theme/color';
import { Point } from '../../../geometry/point/point';
import { PointBuilder } from '../../../geometry/point/point-builder';

export const DefaultGaussianDensityPointProps = {
    ...DefaultUnitsPointProps,

    sizeTheme: { name: 'uniform', value: 1.5 } as SizeThemeProps,
    colorTheme: { name: 'uniform', value: Color(0xAAAAAA) } as ColorThemeProps,
    pointSizeAttenuation: false,

    resolutionFactor: 7,
    probeRadius: 0,
    isoValue: 1.5,
}
export type GaussianDensityPointProps = typeof DefaultGaussianDensityPointProps

export interface GaussianDensityProps {
    sizeTheme: SizeThemeProps

    resolutionFactor: number
    probeRadius: number
    isoValue: number
}

export async function createGaussianDensityPoint(ctx: RuntimeContext, unit: Unit, structure: Structure, props: GaussianDensityProps, point?: Point) {
    const { transform, field: { space, data } } = await computeGaussianDensity(unit, structure, props).runAsChild(ctx)

    const { dimensions, get } = space
    const [ xn, yn, zn ] = dimensions

    const n = xn * yn * zn * 3
    const builder = PointBuilder.create(n, n / 10, point)

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
                if (i % 100000 === 0 && ctx.shouldUpdate) {
                    await ctx.update({ message: 'Creating density points', current: i, max: n });
                }
                ++i
            }
        }
    }
    return builder.getPoint()
}

export function GaussianDensityPointVisual(): UnitsVisual<GaussianDensityPointProps> {
    return UnitsPointVisual<GaussianDensityPointProps>({
        defaultProps: DefaultGaussianDensityPointProps,
        createPoint: createGaussianDensityPoint,
        createLocationIterator: StructureElementIterator.fromGroup,
        getLoci: () => EmptyLoci,
        mark: () => false,
        setUpdateState: (state: VisualUpdateState, newProps: GaussianDensityPointProps, currentProps: GaussianDensityPointProps) => {

        }
    })
}