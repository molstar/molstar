/**
 * Copyright (c) 2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 *
 * adapted from https://github.com/internalfx/distinct-colors (ISC License Copyright (c) 2015, InternalFX Inc.)
 * which is heavily inspired by http://tools.medialab.sciences-po.fr/iwanthue/
 */

import { Lab } from './spaces/lab';
import { Hcl } from './spaces/hcl';
import { deepClone } from '../../mol-util/object';
import { deepEqual } from '../../mol-util';
import { arraySum } from '../../mol-util/array';
import { ParamDefinition as PD } from '../../mol-util/param-definition';

export const DistinctColorsParams = {
    hue: PD.Interval([1, 360], { min: 0, max: 360, step: 1 }),
    chroma: PD.Interval([40, 70], { min: 0, max: 100, step: 1 }),
    luminance: PD.Interval([15, 85], { min: 0, max: 100, step: 1 }),

    clusteringStepCount: PD.Numeric(50, { min: 10, max: 200, step: 1 }, { isHidden: true }),
    minSampleCount: PD.Numeric(800, { min: 100, max: 5000, step: 100 }, { isHidden: true })
};
export type DistinctColorsParams = typeof DistinctColorsParams
export type DistinctColorsProps = PD.Values<typeof DistinctColorsParams>

function distance(colorA: Lab, colorB: Lab) {
    return Math.sqrt(
        Math.pow(Math.abs(colorA[0] - colorB[0]), 2) +
        Math.pow(Math.abs(colorA[1] - colorB[1]), 2) +
        Math.pow(Math.abs(colorA[2] - colorB[2]), 2)
    );
}

const LabTolerance = 2;
const tmpCheckColorHcl = [0, 0, 0] as Hcl;
const tmpCheckColorLab = [0, 0, 0] as Lab;
function checkColor(lab: Lab, props: DistinctColorsProps) {
    Lab.toHcl(tmpCheckColorHcl, lab);
    // roundtrip to RGB for conversion tolerance testing
    Lab.fromColor(tmpCheckColorLab, Lab.toColor(lab));

    return (
        tmpCheckColorHcl[0] >= props.hue[0] &&
        tmpCheckColorHcl[0] <= props.hue[1] &&
        tmpCheckColorHcl[1] >= props.chroma[0] &&
        tmpCheckColorHcl[1] <= props.chroma[1] &&
        tmpCheckColorHcl[2] >= props.luminance[0] &&
        tmpCheckColorHcl[2] <= props.luminance[1] &&
        tmpCheckColorLab[0] >= (lab[0] - LabTolerance) &&
        tmpCheckColorLab[0] <= (lab[0] + LabTolerance) &&
        tmpCheckColorLab[1] >= (lab[1] - LabTolerance) &&
        tmpCheckColorLab[1] <= (lab[1] + LabTolerance) &&
        tmpCheckColorLab[2] >= (lab[2] - LabTolerance) &&
        tmpCheckColorLab[2] <= (lab[2] + LabTolerance)
    );
}

function sortByContrast(colors: Lab[]) {
    const unsortedColors = colors.slice(0);
    const sortedColors = [unsortedColors.shift()!];
    while (unsortedColors.length > 0) {
        const lastColor = sortedColors[sortedColors.length - 1];
        let nearest = 0;
        let maxDist = Number.MIN_SAFE_INTEGER;
        for (let i = 0; i < unsortedColors.length; ++i) {
            const dist = distance(lastColor, unsortedColors[i]);
            if (dist > maxDist) {
                maxDist = dist;
                nearest = i;
            }
        }
        sortedColors.push(unsortedColors.splice(nearest, 1)[0]);
    }
    return sortedColors;
}

function getSamples(count: number, p: DistinctColorsProps) {
    const samples = new Map<string, Lab>();
    const rangeDivider = Math.cbrt(count) * 1.001;

    const hStep = (p.hue[1] - p.hue[0]) / rangeDivider;
    const cStep = (p.chroma[1] - p.chroma[0]) / rangeDivider;
    const lStep = (p.luminance[1] - p.luminance[0]) / rangeDivider;
    for (let h = p.hue[0]; h <= p.hue[1]; h += hStep) {
        for (let c = p.chroma[0]; c <= p.chroma[1]; c += cStep) {
            for (let l = p.luminance[0]; l <= p.luminance[1]; l += lStep) {
                const lab = Lab.fromHcl(Lab(), Hcl.create(h, c, l));
                if (checkColor(lab, p)) samples.set(lab.toString(), lab);
            }
        }
    }

    return Array.from(samples.values());
}

/**
 * Create a list of visually distinct colors
 */
export function distinctColors(count: number, props: Partial<DistinctColorsProps> = {}) {
    const p = { ...PD.getDefaultValues(DistinctColorsParams), ...props };

    if (count <= 0) return [];

    const samples = getSamples(Math.max(p.minSampleCount, count * 5), p);
    if (samples.length < count) {
        throw new Error('Not enough samples to generate distinct colors, increase sample count.');
    }

    const colors: Lab[] = [];
    const zonesProto: (Lab[])[] = [];
    const sliceSize = Math.floor(samples.length / count);

    for (let i = 0; i < samples.length; i += sliceSize) {
        colors.push(samples[i]);
        zonesProto.push([]);
        if (colors.length >= count) break;
    }

    for (let step = 1; step <= p.clusteringStepCount; ++step) {
        const zones = deepClone(zonesProto);

        // Find closest color for each sample
        for (let i = 0; i < samples.length; ++i) {
            let minDist = Number.MAX_SAFE_INTEGER;
            let nearest = 0;
            for (let j = 0; j < colors.length; j++) {
                const dist = distance(samples[i], colors[j]);
                if (dist < minDist) {
                    minDist = dist;
                    nearest = j;
                }
            }
            zones[nearest].push(samples[i]);
        }

        const lastColors = deepClone(colors);

        for (let i = 0; i < zones.length; ++i) {
            const zone = zones[i];
            const size = zone.length;
            const Ls: number[] = [];
            const As: number[] = [];
            const Bs: number[] = [];

            for (let sample of zone) {
                Ls.push(sample[0]);
                As.push(sample[1]);
                Bs.push(sample[2]);
            }

            const lAvg = arraySum(Ls) / size;
            const aAvg = arraySum(As) / size;
            const bAvg = arraySum(Bs) / size;

            colors[i] = [lAvg, aAvg, bAvg] as Lab;
        }

        if (deepEqual(lastColors, colors)) break;
    }

    return sortByContrast(colors).map(c => Lab.toColor(c));
}