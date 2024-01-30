/**
 * Copyright (c) 2019-24 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Loci } from '../../../mol-model/loci';
import { RuntimeContext } from '../../../mol-task';
import { Text } from '../../../mol-geo/geometry/text/text';
import { ParamDefinition as PD } from '../../../mol-util/param-definition';
import { ShapeRepresentation } from '../representation';
import { Representation, RepresentationParamsGetter, RepresentationContext } from '../../representation';
import { Shape } from '../../../mol-model/shape';
import { TextBuilder } from '../../../mol-geo/geometry/text/text-builder';
import { Sphere3D } from '../../../mol-math/geometry';
import { lociLabel } from '../../../mol-theme/label';
import { LociLabelTextParams } from './common';

export interface LabelData {
    infos: { loci: Loci, label?: string }[]
}

const TextParams = {
    ...LociLabelTextParams,
};
type TextParams = typeof TextParams

const LabelVisuals = {
    'text': (ctx: RepresentationContext, getParams: RepresentationParamsGetter<LabelData, TextParams>) => ShapeRepresentation(getTextShape, Text.Utils),
};

export const LabelParams = {
    ...TextParams,
    scaleByRadius: PD.Boolean(true),
    visuals: PD.MultiSelect(['text'], PD.objectToOptions(LabelVisuals)),
    snapshotKey: PD.Text('', { isEssential: true, disableInteractiveUpdates: true, description: 'Activate the snapshot with the provided key when clicking on the label' }),
    tooltip: PD.Text('', { isEssential: true, multiline: true, disableInteractiveUpdates: true, placeholder: 'Tooltip', description: 'Tooltip text to be displayed when hovering over the label' }),
};

export type LabelParams = typeof LabelParams
export type LabelProps = PD.Values<LabelParams>

//

const tmpSphere = Sphere3D();

function label(info: { loci: Loci, label?: string }, condensed = false) {
    return info.label || lociLabel(info.loci, { hidePrefix: true, htmlStyling: false, condensed });
}

function getLabelName(data: LabelData) {
    return data.infos.length === 1 ? label(data.infos[0]) : `${data.infos.length} Labels`;
}

//

function buildText(data: LabelData, props: LabelProps, text?: Text): Text {
    const builder = TextBuilder.create(props, 128, 64, text);
    const customLabel = props.customText.trim();
    for (let i = 0, il = data.infos.length; i < il; ++i) {
        const info = data.infos[i];
        const sphere = Loci.getBoundingSphere(info.loci, tmpSphere);
        if (!sphere) continue;
        const { center, radius } = sphere;
        const text = customLabel || label(info, true);
        builder.add(text, center[0], center[1], center[2], props.scaleByRadius ? radius / 0.9 : 0, props.scaleByRadius ? Math.max(1, radius) : 1, i);
    }
    return builder.getText();
}

function getTextShape(ctx: RuntimeContext, data: LabelData, props: LabelProps, shape?: Shape<Text>) {
    const text = buildText(data, props, shape && shape.geometry);
    const name = getLabelName(data);
    const tooltip = props.tooltip?.trim() ?? '';
    const customLabel = props.customText.trim();
    let getLabel: (groupId: number) => any;

    if (tooltip) {
        getLabel = (_: number) => tooltip;
    } else if (customLabel) {
        getLabel = (_: number) => customLabel;
    } else {
        getLabel = (groupId: number) => label(data.infos[groupId]);
    }

    return Shape.create(name, data, text, () => props.textColor, () => props.textSize, getLabel);
}

//

export type LabelRepresentation = Representation<LabelData, LabelParams>
export function LabelRepresentation(ctx: RepresentationContext, getParams: RepresentationParamsGetter<LabelData, LabelParams>): LabelRepresentation {
    return Representation.createMulti('Label', ctx, getParams, Representation.StateBuilder, LabelVisuals as unknown as Representation.Def<LabelData, LabelParams>);
}