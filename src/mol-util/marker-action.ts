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

    export function isReverse(a: MarkerAction, b: MarkerAction) {
        return (
            (a === MarkerAction.Highlight && b === MarkerAction.RemoveHighlight) ||
            (a === MarkerAction.RemoveHighlight && b === MarkerAction.Highlight) ||
            (a === MarkerAction.Select && b === MarkerAction.Deselect) ||
            (a === MarkerAction.Deselect && b === MarkerAction.Select) ||
            (a === MarkerAction.Toggle && b === MarkerAction.Toggle)
        );
    }
}

export function setMarkerValue(array: Uint8Array, status: 0 | 1 | 2 | 3, count: number) {
    array.fill(status, 0, count);
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
        const viewStart = (start + 3) >> 2;
        const viewEnd = viewStart + ((end - 4 * viewStart) >> 2);

        if (viewEnd <= viewStart) {
            // avoid edge cases with overlapping front/end intervals
            for (let i = start; i < end; ++i) {
                applyMarkerActionAtPosition(array, i, action);
            }
            return true;
        }

        const view = new Uint32Array(array.buffer, 0, array.buffer.byteLength >> 2);

        const frontStart = start;
        const frontEnd = Math.min(4 * viewStart, end);
        const backStart = Math.max(start, 4 * viewEnd);
        const backEnd = end;

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
                for (let i = viewStart; i < viewEnd; ++i) view[i] = 0;
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


export interface MarkerInfo {
    /**
     * 0: none marked;
     * 1: all marked;
     * -1: unclear, need to be calculated
     */
    average: 0 | 1 | -1
    /**
     * 0: none marked;
     * 1: all highlighted;
     * 2: all selected;
     * 3: all highlighted and selected
     * -1: mixed/unclear
     */
    status: 0 | 1 | 2 | 3 | -1
}

export function getMarkerInfo(action: MarkerAction, currentStatus: MarkerInfo['status']): MarkerInfo {
    let average: MarkerInfo['average'] = -1;
    let status: MarkerInfo['status'] = -1;
    switch (action) {
        case MarkerAction.Highlight:
            if (currentStatus === 0 || currentStatus === 1) {
                average = 1;
                status = 1;
            } else if (currentStatus === 2 || currentStatus === 3) {
                average = 1;
                status = 3;
            } else {
                average = 1;
            }
            break;
        case MarkerAction.RemoveHighlight:
            if (currentStatus === 0 || currentStatus === 1) {
                average = 0;
                status = 0;
            } else if (currentStatus === 2 || currentStatus === 3) {
                average = 1;
                status = 2;
            }
            break;
        case MarkerAction.Select:
            if (currentStatus === 1 || currentStatus === 3) {
                average = 1;
                status = 3;
            } else if (currentStatus === 0 || currentStatus === 2) {
                average = 1;
                status = 2;
            } else {
                average = 1;
            }
            break;
        case MarkerAction.Deselect:
            if (currentStatus === 1 || currentStatus === 3) {
                average = 1;
                status = 1;
            } else if (currentStatus === 0 || currentStatus === 2) {
                average = 0;
                status = 0;
            }
            break;
        case MarkerAction.Toggle:
            if (currentStatus === 1) {
                average = 1;
                status = 3;
            } else if (currentStatus === 2) {
                average = 0;
                status = 0;
            } else if (currentStatus === 3) {
                average = 1;
                status = 1;
            } else if (currentStatus === 0) {
                average = 1;
                status = 2;
            }
            break;
        case MarkerAction.Clear:
            average = 0;
            status = 0;
            break;
    }
    return { average, status };
}

/**
 * Assumes the action is applied to a partial set that is
 * neither the empty set nor the full set.
 */
export function getPartialMarkerAverage(action: MarkerAction, currentStatus: MarkerInfo['status']) {
    switch (action) {
        case MarkerAction.Highlight:
            return 0.5;
        case MarkerAction.RemoveHighlight:
            if (currentStatus === 0) {
                return 0;
            } else if (currentStatus === 2 || currentStatus === 3) {
                return 0.5;
            } else { // 1 | -1
                return -1;
            }
        case MarkerAction.Select:
            return 0.5;
        case MarkerAction.Deselect:
            if (currentStatus === 1 || currentStatus === 3) {
                return 0.5;
            } else if (currentStatus === 0) {
                return 0;
            } else { // 2 | -1
                return -1;
            }
        case MarkerAction.Toggle:
            if (currentStatus === -1) {
                return -1;
            } else { // 0 | 1 | 2 | 3
                return 0.5;
            }
        case MarkerAction.Clear:
            if (currentStatus === -1) {
                return -1;
            } else if (currentStatus === 0) {
                return 0;
            } else { // 1 | 2 | 3
                return 0.5;
            }
        case MarkerAction.None:
            return -1;
        default:
            assertUnreachable(action);
    }
}