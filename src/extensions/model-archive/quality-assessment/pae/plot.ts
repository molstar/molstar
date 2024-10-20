/**
 * Copyright (c) 2024 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Model, ResidueIndex } from '../../../../mol-model/structure';
import { AtomicHierarchy } from '../../../../mol-model/structure/model/properties/atomic';
import { Color } from '../../../../mol-util/color';
import { QualityAssessment, QualityAssessmentProvider } from '../prop';


const DefaultPAEColorRange = [0x00441B, 0xF7FCF5] as [Color, Color];

type ResidueRangeInfo = { startOffset: number, endOffset: number, label: string };

function drawAssessmentPNG(model: Model, assessment: QualityAssessment, name: string, colorRange: [Color, Color], noDataColor: Color) {
    const metric = assessment.localPairwiseMetrics.get(name)!;
    const info = assessment.localPairwiseMetricInfo.get(name)!;

    const [minColor, maxColor] = colorRange;
    const range = info.maxResidueIndex - info.minResidueIndex;
    const valueRange = info.maxMetric - info.minMetric;

    const canvas = document.createElement('canvas');
    canvas.width = range;
    canvas.height = range;
    const ctx = canvas.getContext('2d')!;
    ctx.fillStyle = Color.toStyle(noDataColor);
    ctx.fillRect(0, 0, canvas.width, canvas.height);

    metric.forEach((other, a) => {
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

    const chains: ResidueRangeInfo[] = [];
    const hierarchy = model.atomicHierarchy;
    const { label_asym_id } = hierarchy.chains;

    let cI = AtomicHierarchy.residueChainIndex(hierarchy, info.minResidueIndex as ResidueIndex);
    let currentChain: ResidueRangeInfo = { startOffset: 0, endOffset: 1, label: label_asym_id.value(cI) };
    chains.push(currentChain);

    for (let i = 1; i < range; i++) {
        cI = AtomicHierarchy.residueChainIndex(hierarchy, (info.minResidueIndex + i) as ResidueIndex);
        const asym_id = label_asym_id.value(cI);
        if (asym_id === currentChain.label) {
            currentChain.endOffset = i + 1;
        } else {
            currentChain = { startOffset: i, endOffset: i + 1, label: asym_id };
            chains.push(currentChain);
        }
    }

    return {
        model,
        metric,
        info,
        chains,
        colorRange: [Color.toStyle(colorRange[0]), Color.toStyle(colorRange[1])] as const,
        png: canvas.toDataURL('png')
    };
}

export function drawPAEPng(model: Model) {
    const assessment = QualityAssessmentProvider.get(model).value;
    if (!assessment || !assessment.localPairwiseMetrics.has('PAE')) return undefined;
    return drawAssessmentPNG(model, assessment, 'PAE', DefaultPAEColorRange, Color(0xdddddd));
}

export type PAEDrawing = ReturnType<typeof drawPAEPng>