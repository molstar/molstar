/**
 * Copyright (c) 2018-2025 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author Lukáš Polák <admin@lukaspolak.cz>
 */

import { Color, ColorListEntry } from './color';
import { getColorListFromName, ColorListName } from './lists';
import { defaults } from '../../mol-util';
import { NumberArray } from '../../mol-util/type-helpers';
import { ScaleLegend } from '../legend';
import { SortedArray } from '../../mol-data/int';
import { clamp } from '../../mol-math/interpolate';

export interface ColorScale {
    /** Returns hex color for given value */
    color: (value: number) => Color
    /** Copies color to rgb int8 array */
    colorToArray: (value: number, array: NumberArray, offset: number) => void
    /** Copies normalized (0 to 1) hex color to rgb array */
    normalizedColorToArray: (value: number, array: NumberArray, offset: number) => void
    /**  */
    setDomain: (min: number, max: number) => void
    /** Legend */
    readonly legend: ScaleLegend
}

export const DefaultColorScaleProps = {
    domain: [0, 1] as [number, number],
    reverse: false,
    listOrName: 'red-yellow-blue' as ColorListEntry[] | ColorListName,
    minLabel: '' as string | undefined,
    maxLabel: '' as string | undefined,
};
export type ColorScaleProps = Partial<typeof DefaultColorScaleProps>

type ColorScaleType = 'continuous' | 'discrete'

export namespace ColorScale {
    export function create(props: ColorScaleProps): ColorScale {
        return createColorScaleByType(props, 'continuous');
    }

    export function createDiscrete(props: ColorScaleProps): ColorScale {
        return createColorScaleByType(props, 'discrete');
    }

    function createColorScaleByType(props: ColorScaleProps, type: ColorScaleType): ColorScale {
        const { domain, reverse, listOrName } = { ...DefaultColorScaleProps, ...props };
        const list = typeof listOrName === 'string' ? getColorListFromName(listOrName).list : listOrName;

        const colors = reverse ? list.slice().reverse() : list;

        let diff = 0, min = 0, max = 0;
        function setDomain(_min: number, _max: number) {
            min = _min;
            max = _max;
            diff = (max - min) || 1;
        }
        setDomain(domain[0], domain[1]);

        const minLabel = defaults(props.minLabel, min.toString());
        const maxLabel = defaults(props.maxLabel, max.toString());

        let color: (v: number) => Color;

        const hasOffsets = colors.every(c => Array.isArray(c));
        if (hasOffsets) {
            const sorted = [...colors] as [Color, number][];
            sorted.sort((a, b) => a[1] - b[1]);

            const src = sorted.map(c => c[0]);
            const off = SortedArray.ofSortedArray(sorted.map(c => c[1]));
            const max = src.length - 1;

            switch (type) {
                case 'continuous': color = (value: number) => valueToColorWithOffsets(value, src, off, min, max, diff); break;
                case 'discrete': color = (value: number) => valueToDiscreteColorWithOffsets(value, src, off, min, max, diff); break;
            }
        } else {
            switch (type) {
                case 'continuous': color = (value: number) => valueToColor(value, colors, min, max, diff); break;
                case 'discrete': color = (value: number) => valueToDiscreteColor(value, colors, min, max, diff); break;
            }
        }
        return {
            color,
            colorToArray: (value: number, array: NumberArray, offset: number) => {
                Color.toArray(color(value), array, offset);
            },
            normalizedColorToArray: (value: number, array: NumberArray, offset: number) => {
                Color.toArrayNormalized(color(value), array, offset);
            },
            setDomain,
            get legend() { return ScaleLegend(minLabel, maxLabel, colors); }
        };
    }

    function valueToColorWithOffsets(value: number, src: Color[], off: SortedArray<number>, min: number, max: number, diff: number) {
        const t = clamp((value - min) / diff, 0, 1);
        const i = SortedArray.findPredecessorIndex(off, t);

        if (i === 0) {
            return src[min];
        } else if (i > max) {
            return src[max];
        }

        const o1 = off[i - 1], o2 = off[i];
        const t1 = clamp((t - o1) / (o2 - o1), 0, 1); // TODO: cache the deltas?

        return Color.interpolate(src[i - 1], src[i], t1);
    }

    function valueToColor(value: number, colors: ColorListEntry[], min: number, max: number, diff: number) {
        const t = Math.min(colors.length - 1, Math.max(0, ((value - min) / diff) * colors.length - 1));
        const tf = Math.floor(t);
        const c1 = colors[tf] as Color;
        const c2 = colors[Math.ceil(t)] as Color;
        return Color.interpolate(c1, c2, t - tf);
    }

    function valueToDiscreteColorWithOffsets(value: number, src: Color[], off: SortedArray<number>, min: number, max: number, diff: number) {
        if (src.length === 0) {
            return Color.fromRgb(0, 0, 0);
        }

        const t = clamp((value - min) / diff, 0, 1);
        const i = SortedArray.findPredecessorIndex(off, t);

        if (i === 0) {
            return src[min] as Color;
        } else if (i > max) {
            return src[max] as Color;
        }

        return src[i] as Color;
    }

    function valueToDiscreteColor(value: number, colors: ColorListEntry[], min: number, max: number, diff: number) {
        if (colors.length === 0) {
            return Color.fromRgb(0, 0, 0);
        }

        const intervalSize = diff / colors.length;

        if (value <= min) {
            return colors[0] as Color;
        } else if (value >= max) {
            return colors[colors.length - 1] as Color;
        }

        const i = Math.min(colors.length - 1, Math.floor((value - min) / intervalSize));

        return colors[i] as Color;
    }
}
