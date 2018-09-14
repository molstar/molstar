/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Unit, Structure } from 'mol-model/structure';
import { RuntimeContext } from 'mol-task'
import { UnitsVisual, VisualUpdateState } from '..';
import { getElementLoci, StructureElementIterator, markElement } from './util/element';
import { Vec3 } from 'mol-math/linear-algebra';
import { SizeThemeProps } from 'mol-view/theme/size';
import { UnitsPointVisual, DefaultUnitsPointProps } from '../units-visual';
import { Point } from '../../../geometry/point/point';
import { PointBuilder } from '../../../geometry/point/point-builder';

export const DefaultElementPointProps = {
    ...DefaultUnitsPointProps,

    sizeTheme: { name: 'uniform', value: 0.2 } as SizeThemeProps,
    pointSizeAttenuation: true,
}
export type ElementPointProps = typeof DefaultElementPointProps

// TODO size

export async function createElementPoint(ctx: RuntimeContext, unit: Unit, structure: Structure, props: ElementPointProps, point: Point) {
    const elements = unit.elements
    const n = elements.length
    const builder = PointBuilder.create(n, n / 10, point)

    const pos = unit.conformation.invariantPosition
    const p = Vec3.zero()

    for (let i = 0; i < n; ++i) {
        pos(elements[i], p)
        builder.add(p[0], p[1], p[2], i)

        if (i % 10000 === 0 && ctx.shouldUpdate) {
            await ctx.update({ message: 'Creating points', current: i, max: n });
        }
    }
    return builder.getPoint()
}

export function ElementPointVisual(): UnitsVisual<ElementPointProps> {
    return UnitsPointVisual<ElementPointProps>({
        defaultProps: DefaultElementPointProps,
        createPoint: createElementPoint,
        createLocationIterator: StructureElementIterator.fromGroup,
        getLoci: getElementLoci,
        mark: markElement,
        setUpdateState: (state: VisualUpdateState, newProps: ElementPointProps, currentProps: ElementPointProps) => {

        }
    })
}