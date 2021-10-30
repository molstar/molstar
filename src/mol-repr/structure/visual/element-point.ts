/**
 * Copyright (c) 2018-2021 mol* contributors, licensed under MIT, See LICENSE file for more info.
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
import { ElementIterator, getElementLoci, eachElement, makeElementIgnoreTest } from './util/element';
import { VisualUpdateState } from '../../util';
import { Sphere3D } from '../../../mol-math/geometry';

export const ElementPointParams = {
    ...UnitsPointsParams,
    pointSizeAttenuation: PD.Boolean(false),
    ignoreHydrogens: PD.Boolean(false),
    traceOnly: PD.Boolean(false),
};
export type ElementPointParams = typeof ElementPointParams

// TODO size

export function createElementPoint(ctx: VisualContext, unit: Unit, structure: Structure, theme: Theme, props: PD.Values<ElementPointParams>, points: Points) {
    // TODO sizeFactor

    const { child } = structure;
    if (child && !child.unitMap.get(unit.id)) return Points.createEmpty(points);

    const elements = unit.elements;
    const n = elements.length;
    const builder = PointsBuilder.create(n, n / 10, points);

    const p = Vec3();
    const pos = unit.conformation.invariantPosition;
    const ignore = makeElementIgnoreTest(structure, unit, props);

    if (ignore) {
        for (let i = 0; i < n; ++i) {
            if (ignore(elements[i])) continue;
            pos(elements[i], p);
            builder.add(p[0], p[1], p[2], i);
        }
    } else {
        for (let i = 0; i < n; ++i) {
            pos(elements[i], p);
            builder.add(p[0], p[1], p[2], i);
        }
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
            state.createGeometry = (
                newProps.ignoreHydrogens !== currentProps.ignoreHydrogens ||
                newProps.traceOnly !== currentProps.traceOnly
            );
        }
    }, materialId);
}