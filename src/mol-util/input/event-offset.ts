/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

/*
 * This code has been modified from https://github.com/mattdesl/mouse-event-offset,
 * copyright (c) 2014 Matt DesLauriers. MIT License
 */

import { Vec2 } from '../../mol-math/linear-algebra';

const rootPosition = { left: 0, top: 0 };

export function eventOffset (out: Vec2, ev: MouseEvent | Touch, target: Element) {
    const cx = ev.clientX || 0;
    const cy = ev.clientY || 0;
    const rect = getBoundingClientOffset(target);
    out[0] = cx - rect.left;
    out[1] = cy - rect.top;
    return out;
}

function getBoundingClientOffset (element: Element | Window | Document) {
    if (element instanceof Window || element instanceof Document || element === document.body) {
        return rootPosition;
    } else {
        return (element as Element).getBoundingClientRect();
    }
}