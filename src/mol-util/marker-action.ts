/**
 * Copyright (c) 2018-2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { OrderedSet, Interval } from '../mol-data/int';
import BitFlags from './bit-flags';

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
    let v = array[i];
    switch (action) {
        case MarkerAction.Highlight:
            if (v % 2 === 0) {
                array[i] = v + 1;
                return true;
            }
            return false;
        case MarkerAction.RemoveHighlight:
            if (v % 2 !== 0) {
                array[i] = v - 1;
                return true;
            }
            return false;
        case MarkerAction.Select:
            if (v < 2) {
                array[i] = v + 2;
                return true;
            }
            return false;
        case MarkerAction.Deselect:
            array[i] = v % 2;
            return array[i] !== v;
        case MarkerAction.Toggle:
            if (v >= 2) array[i] = v - 2;
            else array[i] = v + 2;
            return true;
        case MarkerAction.Clear:
            array[i] = 0;
            return v !== 0;
    }
    return false;
}

export function applyMarkerAction(array: Uint8Array, set: OrderedSet, action: MarkerAction) {
    let changed = false;
    if (Interval.is(set)) {
        for (let i = Interval.start(set), _i = Interval.end(set); i < _i; i++) {
            changed = applyMarkerActionAtPosition(array, i, action) || changed;
        }
    } else {
        for (let i = 0, _i = set.length; i < _i; i++) {
            changed = applyMarkerActionAtPosition(array, set[i], action) || changed;
        }
    }
    return changed;
}
