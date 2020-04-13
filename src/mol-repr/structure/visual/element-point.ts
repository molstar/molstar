/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { ParamDefinition as PD } from '../../../mol-util/param-definition';
import { UnitsPointsParams, UnitsVisual, UnitsPointsVisual } from '../units-visual';
import { VisualContext } from '../../visual';
import { Unit, Structure } from '../../../mol-model/structure';
import { Theme } from '../../../mol-theme/theme';
import { Points } from '../../../mol-geo/geometry/points/points';
import { PointsBuilder } from '../../../mol-geo/geometry/points/points-builder';
import { Vec3 } from '../../../mol-math/linear-algebra';
import { ElementIterator, getElementLoci, eachElement } from './util/element';
import { VisualUpdateState } from '../../util';
import { Sphere3D } from '../../../mol-math/geometry';

export const ElementPointParams = {
    ...UnitsPointsParams,
    // sizeFactor: PD.Numeric(1.0, { min: 0, max: 10, step: 0.01 }),
    pointSizeAttenuation: PD.Boolean(false),
    showHydrogens: PD.Boolean(true),
};
export type ElementPointParams = typeof ElementPointParams

// TODO size

export function createElementPoint(ctx: VisualContext, unit: Unit, structure: Structure, theme: Theme, props: PD.Values<ElementPointParams>, points: Points) {
    // const { sizeFactor } = props

    const elements = unit.elements;
    const n = elements.length;
    const builder = PointsBuilder.create(n, n / 10, points);

    const pos = unit.conformation.invariantPosition;
    const p = Vec3.zero();

    for (let i = 0; i < n; ++i) {
        pos(elements[i], p);
        builder.add(p[0], p[1], p[2], i);
    }

    const pt = builder.getPoints();

    const sphere = Sphere3D.expand(Sphere3D(), unit.boundary.sphere, 1 * props.sizeFactor);
    pt.setBoundingSphere(sphere);

    return pt;
}

export function ElementPointVisual(materialId: number): UnitsVisual<ElementPointParams> {
    return UnitsPointsVisual<ElementPointParams>({
        defaultProps: PD.getDefaultValues(ElementPointParams),
        createGeometry: createElementPoint,
        createLocationIterator: ElementIterator.fromGroup,
        getLoci: getElementLoci,
        eachLocation: eachElement,
        setUpdateState: (state: VisualUpdateState, newProps: PD.Values<ElementPointParams>, currentProps: PD.Values<ElementPointParams>) => {

        }
    }, materialId);
}