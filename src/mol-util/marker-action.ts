/**
 * Copyright (c) 2018-2021 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { OrderedSet, Interval } from '../mol-data/int';
import { BitFlags } from './bit-flags';
import { assertUnreachable } from './type-helpers';

export enum MarkerAction {
    None = 0x0,
    Highlight = 0x1,
    RemoveHighlight = 0x2,
    Select = 0x4,
    Deselect = 0x8,
    Toggle = 0x10,
    Clear = 0x20
}

export type MarkerActions = BitFlags<MarkerAction>
export namespace MarkerActions {
    export const is: (m: MarkerActions, f: MarkerAction) => boolean = BitFlags.has;

    export const All = (
        MarkerAction.Highlight | MarkerAction.RemoveHighlight |
        MarkerAction.Select | MarkerAction.Deselect | MarkerAction.Toggle |
        MarkerAction.Clear
    ) as MarkerActions;
    export const Highlighting = (
        MarkerAction.Highlight | MarkerAction.RemoveHighlight |
        MarkerAction.Clear
    ) as MarkerActions;
    export const Selecting = (
        MarkerAction.Select | MarkerAction.Deselect | MarkerAction.Toggle |
        MarkerAction.Clear
    ) as MarkerActions;
}

export function applyMarkerActionAtPosition(array: Uint8Array, i: number, action: MarkerAction) {
    switch (action) {
        case MarkerAction.Highlight: array[i] |= 1; break;
        case MarkerAction.RemoveHighlight: array[i] &= ~1; break;
        case MarkerAction.Select: array[i] |= 2; break;
        case MarkerAction.Deselect: array[i] &= ~2; break;
        case MarkerAction.Toggle: array[i] ^= 2; break;
        case MarkerAction.Clear: array[i] = 0; break;
    }
}

export function applyMarkerAction(array: Uint8Array, set: OrderedSet, action: MarkerAction) {
    if (action === MarkerAction.None) return false;

    if (Interval.is(set)) {
        const start = Interval.start(set);
        const end = Interval.end(set);
        const view = new Uint32Array(array.buffer, 0, Math.floor(array.buffer.byteLength / 4));
        const viewStart = Math.ceil(start / 4);
        const viewEnd = Math.min(view.length, Math.floor(end / 4));

        const middleStart = viewStart * 4;
        const middleEnd = viewEnd * 4;
        const frontStart = start;
        const frontEnd = frontStart === middleStart ? frontStart : middleStart;
        const backEnd = end;
        const backStart = backEnd === middleEnd ? backEnd : middleEnd;

        switch (action) {
            case MarkerAction.Highlight:
                for (let i = viewStart; i < viewEnd; ++i) view[i] |= 0x01010101;
                break;
            case MarkerAction.RemoveHighlight:
                for (let i = viewStart; i < viewEnd; ++i) view[i] &= ~0x01010101;
                break;
            case MarkerAction.Select:
                for (let i = viewStart; i < viewEnd; ++i) view[i] |= 0x02020202;
                break;
            case MarkerAction.Deselect:
                for (let i = viewStart; i < viewEnd; ++i) view[i] &= ~0x02020202;
                break;
            case MarkerAction.Toggle:
                for (let i = viewStart; i < viewEnd; ++i) view[i] ^= 0x02020202;
                break;
            case MarkerAction.Clear:
                view.fill(0);
                break;
            default:
                assertUnreachable(action);
        }

        for (let i = frontStart; i < frontEnd; ++i) {
            applyMarkerActionAtPosition(array, i, action);
        }

        for (let i = backStart; i < backEnd; ++i) {
            applyMarkerActionAtPosition(array, i, action);
        }
    } else {
        switch (action) {
            case MarkerAction.Highlight:
                for (let i = 0, il = set.length; i < il; ++i) array[set[i]] |= 1;
                break;
            case MarkerAction.RemoveHighlight:
                for (let i = 0, il = set.length; i < il; ++i) array[set[i]] &= ~1;
                break;
            case MarkerAction.Select:
                for (let i = 0, il = set.length; i < il; ++i) array[set[i]] |= 2;
                break;
            case MarkerAction.Deselect:
                for (let i = 0, il = set.length; i < il; ++i) array[set[i]] &= ~2;
                break;
            case MarkerAction.Toggle:
                for (let i = 0, il = set.length; i < il; ++i) array[set[i]] ^= 2;
                break;
            case MarkerAction.Clear:
                for (let i = 0, il = set.length; i < il; ++i) array[set[i]] = 0;
                break;
            default:
                assertUnreachable(action);
        }
    }
    return true;
}
