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
import { SizeThemeOptions, SizeThemeName } from 'mol-theme/size';
import { UnitsPointsVisual, UnitsPointsParams } from '../units-visual';
import { Points } from '../../../geometry/points/points';
import { PointsBuilder } from '../../../geometry/points/points-builder';
import { SelectParam, NumberParam, BooleanParam, paramDefaultValues } from 'mol-util/parameter';

export const ElementPointParams = {
    ...UnitsPointsParams,
    sizeTheme: SelectParam<SizeThemeName>('Size Theme', '', 'uniform', SizeThemeOptions),
    sizeValue: NumberParam('Size Value', '', 3, 0, 20, 0.1),
    pointSizeAttenuation: BooleanParam('Point Size Attenuation', '', false),
}
export const DefaultElementPointProps = paramDefaultValues(ElementPointParams)
export type ElementPointProps = typeof DefaultElementPointProps

// TODO size

export async function createElementPoint(ctx: RuntimeContext, unit: Unit, structure: Structure, props: ElementPointProps, points: Points) {
    const elements = unit.elements
    const n = elements.length
    const builder = PointsBuilder.create(n, n / 10, points)

    const pos = unit.conformation.invariantPosition
    const p = Vec3.zero()

    for (let i = 0; i < n; ++i) {
        pos(elements[i], p)
        builder.add(p[0], p[1], p[2], i)

        if (i % 10000 === 0 && ctx.shouldUpdate) {
            await ctx.update({ message: 'Creating points', current: i, max: n });
        }
    }
    return builder.getPoints()
}

export function ElementPointVisual(): UnitsVisual<ElementPointProps> {
    return UnitsPointsVisual<ElementPointProps>({
        defaultProps: DefaultElementPointProps,
        createGeometry: createElementPoint,
        createLocationIterator: StructureElementIterator.fromGroup,
        getLoci: getElementLoci,
        mark: markElement,
        setUpdateState: (state: VisualUpdateState, newProps: ElementPointProps, currentProps: ElementPointProps) => {

        }
    })
}