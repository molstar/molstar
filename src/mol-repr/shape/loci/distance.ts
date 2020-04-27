/**
 * Copyright (c) 2019-2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Loci } from '../../../mol-model/loci';
import { RuntimeContext } from '../../../mol-task';
import { Lines } from '../../../mol-geo/geometry/lines/lines';
import { Text } from '../../../mol-geo/geometry/text/text';
import { ParamDefinition as PD } from '../../../mol-util/param-definition';
import { ColorNames } from '../../../mol-util/color/names';
import { ShapeRepresentation } from '../representation';
import { Representation, RepresentationParamsGetter, RepresentationContext } from '../../representation';
import { Shape } from '../../../mol-model/shape';
import { LinesBuilder } from '../../../mol-geo/geometry/lines/lines-builder';
import { TextBuilder } from '../../../mol-geo/geometry/text/text-builder';
import { Vec3 } from '../../../mol-math/linear-algebra';
import { MarkerActions, MarkerAction } from '../../../mol-util/marker-action';
import { distanceLabel } from '../../../mol-theme/label';
import { MeasurementRepresentationCommonTextParams } from './common';
import { Sphere3D } from '../../../mol-math/geometry';

export interface DistanceData {
    pairs: Loci.Bundle<2>[]
}

const SharedParams = {
    unitLabel: PD.Text('\u212B', { isEssential: true })
};

const LineParams = {
    ...Lines.Params,
    ...SharedParams,
    lineSizeAttenuation: PD.Boolean(true),
    linesColor: PD.Color(ColorNames.lightgreen, { isEssential: true }),
    linesSize: PD.Numeric(0.075, { min: 0.01, max: 5, step: 0.01 }),
    dashLength: PD.Numeric(0.2, { min: 0.01, max: 0.2, step: 0.01 }),
};
type LineParams = typeof LineParams

const TextParams = {
    ...Text.Params,
    ...SharedParams,
    borderWidth: PD.Numeric(0.2, { min: 0, max: 0.5, step: 0.01 }),
    ...MeasurementRepresentationCommonTextParams
};
type TextParams = typeof TextParams

const DistanceVisuals = {
    'lines': (ctx: RepresentationContext, getParams: RepresentationParamsGetter<DistanceData, LineParams>) => ShapeRepresentation(getLinesShape, Lines.Utils, { modifyState: s => ({ ...s, markerActions: MarkerActions.Highlighting }) }),
    'text': (ctx: RepresentationContext, getParams: RepresentationParamsGetter<DistanceData, TextParams>) => ShapeRepresentation(getTextShape, Text.Utils, { modifyState: s => ({ ...s, markerActions: MarkerAction.None }) }),
};

export const DistanceParams = {
    ...LineParams,
    ...TextParams,
    visuals: PD.MultiSelect(['lines', 'text'], PD.objectToOptions(DistanceVisuals)),
};
export type DistanceParams = typeof DistanceParams
export type DistanceProps = PD.Values<DistanceParams>

//

function getDistanceState() {
    return {
        sphereA: Sphere3D(),
        sphereB: Sphere3D(),

        center: Vec3(),
        distance: 0,
    };
}
type DistanceState = ReturnType<typeof getDistanceState>

function setDistanceState(pair: Loci.Bundle<2>, state: DistanceState) {
    const { sphereA, sphereB, center } = state;

    const [lociA, lociB] = pair.loci;
    Loci.getBoundingSphere(lociA, sphereA);
    Loci.getBoundingSphere(lociB, sphereB);

    Vec3.add(center, sphereA.center, sphereB.center);
    Vec3.scale(center, center, 0.5);
    state.distance = Vec3.distance(sphereA.center, sphereB.center);

    return state;
}

const tmpState = getDistanceState();

function getDistanceName(data: DistanceData, unitLabel: string) {
    return data.pairs.length === 1 ? `Distance ${distanceLabel(data.pairs[0], { unitLabel, measureOnly: true })}` : `${data.pairs.length} Distances`;
}

//

function buildLines(data: DistanceData, props: DistanceProps, lines?: Lines): Lines {
    const builder = LinesBuilder.create(128, 64, lines);
    for (let i = 0, il = data.pairs.length; i < il; ++i) {
        setDistanceState(data.pairs[i], tmpState);
        builder.addFixedLengthDashes(tmpState.sphereA.center, tmpState.sphereB.center, props.dashLength, i);
    }
    return builder.getLines();
}

function getLinesShape(ctx: RuntimeContext, data: DistanceData, props: DistanceProps, shape?: Shape<Lines>) {
    const lines = buildLines(data, props, shape && shape.geometry);
    const name = getDistanceName(data, props.unitLabel);
    const getLabel = (groupId: number ) => distanceLabel(data.pairs[groupId], props);
    return Shape.create(name, data, lines, () => props.linesColor, () => props.linesSize, getLabel);
}

//

function buildText(data: DistanceData, props: DistanceProps, text?: Text): Text {
    const builder = TextBuilder.create(props, 128, 64, text);
    for (let i = 0, il = data.pairs.length; i < il; ++i) {
        setDistanceState(data.pairs[i], tmpState);
        const { center, distance, sphereA, sphereB } = tmpState;
        const label = props.customText || `${distance.toFixed(2)} ${props.unitLabel}`;
        const radius = Math.max(2, sphereA.radius, sphereB.radius);
        const scale = radius / 2;
        builder.add(label, center[0], center[1], center[2], 1, scale, i);
    }
    return builder.getText();
}

function getTextShape(ctx: RuntimeContext, data: DistanceData, props: DistanceProps, shape?: Shape<Text>) {
    const text = buildText(data, props, shape && shape.geometry);
    const name = getDistanceName(data, props.unitLabel);
    const getLabel = (groupId: number ) => distanceLabel(data.pairs[groupId], props);
    return Shape.create(name, data, text, () => props.textColor, () => props.textSize, getLabel);
}

//

export type DistanceRepresentation = Representation<DistanceData, DistanceParams>
export function DistanceRepresentation(ctx: RepresentationContext, getParams: RepresentationParamsGetter<DistanceData, DistanceParams>): DistanceRepresentation {
    return Representation.createMulti('Distance', ctx, getParams, Representation.StateBuilder, DistanceVisuals as unknown as Representation.Def<DistanceData, DistanceParams>);
}