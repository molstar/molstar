/**
 * Copyright (c) 2024 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Model } from '../../../../mol-model/structure';
import { Color } from '../../../../mol-util/color';
import { QualityAssessment, QualityAssessmentProvider } from '../prop';


const DefaultPAEColorRange = [0xF7FCF5 as Color, 0x00441B as Color] as [Color, Color];

function drawAssessmentPNG(assesment: QualityAssessment, name: string, colorRange: [Color, Color]) {
    const data = assesment.localPairwiseMetrics.get(name)!;
    const info = assesment.localPairwiseMetricInfo.get(name)!;

    const [minColor, maxColor] = colorRange;
    const range = info.maxResidueIndex - info.minResidueIndex;
    const valueRange = info.maxMetric - info.minMetric;

    const canvas = document.createElement('canvas');
    canvas.width = range;
    canvas.height = range;
    const ctx = canvas.getContext('2d')!;
    ctx.fillStyle = 'rgb(255, 255, 255)';
    ctx.fillRect(0, 0, canvas.width, canvas.height);

    data.forEach((other, a) => {
        other.forEach((value, b) => {
            const x = a - info.minResidueIndex;
            const y = b - info.minResidueIndex;
            const t = (value - info.minMetric) / valueRange;

            const color = Color.interpolate(minColor, maxColor, t);
            ctx.fillStyle = Color.toStyle(color);
            ctx.fillRect(x, y, 1, 1);
            ctx.fillRect(y, x, 1, 1);
        });
    });

    return canvas.toDataURL('png');
}

export function drawPAEPng(model: Model) {
    const assessment = QualityAssessmentProvider.get(model).value;
    if (!assessment || !assessment.localPairwiseMetrics.has('PAE')) return undefined;
    return drawAssessmentPNG(assessment, 'PAE', DefaultPAEColorRange);
}