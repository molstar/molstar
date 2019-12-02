/**
 * Copyright (c) 2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Loci } from '../../../mol-model/loci';
import { RuntimeContext } from '../../../mol-task';
import { stringToWords } from '../../../mol-util/string';
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

export interface DistanceData {
    pairs: { lociA: Loci, lociB: Loci }[]
}

const LineParams = {
    ...Lines.Params,
    lineSizeAttenuation: PD.Boolean(true),
    linesColor: PD.Color(ColorNames.darkgreen),
    linesSize: PD.Numeric(0.05, { min: 0.01, max: 5, step: 0.01 }),
}
type LineParams = typeof LineParams

const TextParams = {
    ...Text.Params,
    borderWidth: PD.Numeric(0.25, { min: 0, max: 0.5, step: 0.01 }),
    textColor: PD.Color(ColorNames.black),
    textSize: PD.Numeric(0.8, { min: 0.1, max: 5, step: 0.1 }),
    unitLabel: PD.Text('\u212B')
}
type TextParams = typeof TextParams

const DistanceVisuals = {
    'lines': (ctx: RepresentationContext, getParams: RepresentationParamsGetter<DistanceData, LineParams>) => ShapeRepresentation(getLinesShape, Lines.Utils),
    'text': (ctx: RepresentationContext, getParams: RepresentationParamsGetter<DistanceData, TextParams>) => ShapeRepresentation(getTextShape, Text.Utils),
}
type DistanceVisualName = keyof typeof DistanceVisuals
const DistanceVisualOptions = Object.keys(DistanceVisuals).map(name => [name, stringToWords(name)] as [DistanceVisualName, string])

export const DistanceParams = {
    ...LineParams,
    ...TextParams,
    visuals: PD.MultiSelect<DistanceVisualName>(['lines', 'text'], DistanceVisualOptions),
}
export type DistanceParams = typeof DistanceParams
export type DistanceProps = PD.Values<DistanceParams>

//

const tmpStartV = Vec3()
const tmpEndV = Vec3()
const tmpCenterV = Vec3()

function buildLines(data: DistanceData, props: DistanceProps, lines?: Lines): Lines {
    const builder = LinesBuilder.create(128, 64, lines)
    for (let i = 0, il = data.pairs.length; i < il; ++i) {
        const { lociA, lociB } = data.pairs[i]
        Loci.getCenter(lociA, tmpStartV)
        Loci.getCenter(lociB, tmpEndV)
        builder.addFixedLengthDashes(tmpStartV, tmpEndV, 0.1, i)
    }
    return builder.getLines()
}

function getLinesShape(ctx: RuntimeContext, data: DistanceData, props: DistanceProps, shape?: Shape<Lines>) {
    const lines = buildLines(data, props, shape && shape.geometry);
    const getLabel = function (groupId: number ) {
        return 'Distance Line'
    }
    return Shape.create('Distance Lines', data, lines, () => props.linesColor, () => props.linesSize, getLabel)
}

//

function buildText(data: DistanceData, props: DistanceProps, text?: Text): Text {
    const builder = TextBuilder.create(props, 128, 64, text)
    for (let i = 0, il = data.pairs.length; i < il; ++i) {
        const { lociA, lociB } = data.pairs[i]
        Loci.getCenter(lociA, tmpStartV)
        Loci.getCenter(lociB, tmpEndV)
        Vec3.add(tmpCenterV, tmpStartV, tmpEndV)
        Vec3.scale(tmpCenterV, tmpCenterV, 0.5)
        const dist = Vec3.distance(tmpStartV, tmpEndV).toPrecision(2)
        const label = `${dist} ${props.unitLabel}`
        builder.add(label, tmpCenterV[0], tmpCenterV[1], tmpCenterV[2], 0.1, i)
    }
    return builder.getText()
}

function getTextShape(ctx: RuntimeContext, data: DistanceData, props: DistanceProps, shape?: Shape<Text>) {
    const text = buildText(data, props, shape && shape.geometry);
    const getLabel = function (groupId: number ) {
        return 'Distance Text'
    }
    return Shape.create('Distance Text', data, text, () => props.textColor, () => props.textSize, getLabel)
}

//

export type DistanceRepresentation = Representation<DistanceData, DistanceParams>
export function DistanceRepresentation(ctx: RepresentationContext, getParams: RepresentationParamsGetter<DistanceData, DistanceParams>): DistanceRepresentation {
    return Representation.createMulti('Distance', ctx, getParams, Representation.StateBuilder, DistanceVisuals as unknown as Representation.Def<DistanceData, DistanceParams>)
}