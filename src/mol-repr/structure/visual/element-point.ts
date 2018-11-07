/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Unit, Structure } from 'mol-model/structure';
import { UnitsVisual } from '../index';
import { VisualUpdateState } from '../../util';
import { getElementLoci, StructureElementIterator, markElement } from './util/element';
import { Vec3 } from 'mol-math/linear-algebra';
import { SizeThemeOptions, SizeThemeName } from 'mol-theme/size';
import { UnitsPointsVisual, UnitsPointsParams } from '../units-visual';
import { ParamDefinition as PD } from 'mol-util/param-definition';
import { Points } from 'mol-geo/geometry/points/points';
import { PointsBuilder } from 'mol-geo/geometry/points/points-builder';
import { VisualContext } from 'mol-repr';
import { Theme } from 'mol-geo/geometry/geometry';

export const ElementPointParams = {
    ...UnitsPointsParams,
    sizeTheme: PD.Select<SizeThemeName>('Size Theme', '', 'uniform', SizeThemeOptions),
    sizeValue: PD.Numeric('Size Value', '', 3, 0, 20, 0.1),
    pointSizeAttenuation: PD.Boolean('Point Size Attenuation', '', false),
}
export const DefaultElementPointProps = PD.getDefaultValues(ElementPointParams)
export type ElementPointProps = typeof DefaultElementPointProps

// TODO size

export async function createElementPoint(ctx: VisualContext, unit: Unit, structure: Structure, theme: Theme, props: ElementPointProps, points: Points) {
    const elements = unit.elements
    const n = elements.length
    const builder = PointsBuilder.create(n, n / 10, points)

    const pos = unit.conformation.invariantPosition
    const p = Vec3.zero()

    for (let i = 0; i < n; ++i) {
        pos(elements[i], p)
        builder.add(p[0], p[1], p[2], i)

        if (i % 10000 === 0 && ctx.runtime.shouldUpdate) {
            await ctx.runtime.update({ message: 'Creating points', current: i, max: n });
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