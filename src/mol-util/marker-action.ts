import { OrderedSet } from '../mol-data/int';

/**
 * Copyright (c) 2018-2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

export enum MarkerAction {
    Highlight,
    RemoveHighlight,
    Select,
    Deselect,
    Toggle,
    Clear
}

export function applyMarkerAction(array: Uint8Array, set: OrderedSet, action: MarkerAction) {
    let changed = false;
    OrderedSet.forEach(set, i => {
        let v = array[i];
        switch (action) {
            case MarkerAction.Highlight:
                if (v % 2 === 0) {
                    v += 1;
                }
                break;
            case MarkerAction.RemoveHighlight:
                if (v % 2 !== 0) {
                    v -= 1;
                }
                break;
            case MarkerAction.Select:
                if (v < 2)
                    v += 2;
                // v += 2
                break;
            case MarkerAction.Deselect:
                // if (v >= 2) {
                //     v -= 2
                // }
                v = v % 2;
                break;
            case MarkerAction.Toggle:
                if (v >= 2) {
                    v -= 2;
                }
                else {
                    v += 2;
                }
                break;
            case MarkerAction.Clear:
                v = 0;
                break;
        }
        changed = array[i] !== v || changed;
        array[i] = v;
    })
    return changed;
}
