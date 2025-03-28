/**
 * Copyright (c) 2024-25 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Model, ResidueIndex } from '../../../../mol-model/structure';
import { AtomicHierarchy } from '../../../../mol-model/structure/model/properties/atomic';
import { Color } from '../../../../mol-util/color';
import { QualityAssessment } from '../prop';


const DefaultMetricColorRange = [0x00441B, 0xF7FCF5] as [Color, Color];

export type MAResidueRangeInfo = { startOffset: number, endOffset: number, label: string };

function drawMetricPNG(model: Model, metric: QualityAssessment.Pairwise, colorRange: [Color, Color], noDataColor: Color) {
    const [minResidueIndex, maxResidueIndex] = metric.residueRange;
    const [minMetric, maxMetric] = metric.valueRange;
    const [minColor, maxColor] = colorRange;
    const range = maxResidueIndex - minResidueIndex + 1;
    const valueRange = maxMetric - minMetric;
    const values = metric.values;

    const canvas = document.createElement('canvas');
    canvas.width = range;
    canvas.height = range;
    const ctx = canvas.getContext('2d')!;
    ctx.fillStyle = Color.toStyle(noDataColor);
    ctx.fillRect(0, 0, canvas.width, canvas.height);

    const colorCache = new Map<number, string>();
    const getColor = (t: number) => {
        const rounded = Math.round(t * 0xffff);
        if (colorCache.has(rounded)) {
            return colorCache.get(rounded)!;
        }
        const color = Color.interpolate(minColor, maxColor, rounded / 0xffff);
        const style = Color.toStyle(color);
        colorCache.set(rounded, style);
        return style;
    };

    for (let rA = minResidueIndex; rA <= maxResidueIndex; rA++) {
        const row = values[rA];
        if (!row) continue;

        for (let rB = minResidueIndex; rB <= maxResidueIndex; rB++) {
            const value = row[rB];
            if (typeof value !== 'number') continue;

            const x = rA - minResidueIndex;
            const y = rB - minResidueIndex;
            const t = (value - minMetric) / valueRange;
            ctx.fillStyle = getColor(t);
            ctx.fillRect(x, y, 1, 1);

            if (typeof values[rB]?.[rA] !== 'number') {
                ctx.fillRect(y, x, 1, 1);
            }
        }
    }

    const chains: MAResidueRangeInfo[] = [];
    const hierarchy = model.atomicHierarchy;
    const { label_asym_id } = hierarchy.chains;

    let cI = AtomicHierarchy.residueChainIndex(hierarchy, minResidueIndex as ResidueIndex);
    let currentChain: MAResidueRangeInfo = { startOffset: 0, endOffset: 1, label: label_asym_id.value(cI) };
    chains.push(currentChain);

    for (let i = 1; i < range; i++) {
        cI = AtomicHierarchy.residueChainIndex(hierarchy, (minResidueIndex + i) as ResidueIndex);
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
        chains,
        colorRange: [Color.toStyle(colorRange[0]), Color.toStyle(colorRange[1])] as const,
        png: canvas.toDataURL('png')
    };
}

export function maDrawPairwiseMetricPNG(model: Model, metric: QualityAssessment.Pairwise) {
    return drawMetricPNG(model, metric, DefaultMetricColorRange, Color(0xE2E2E2));
}

export type MAPairwiseMetricDrawing = ReturnType<typeof drawMetricPNG>