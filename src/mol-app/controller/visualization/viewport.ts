/*
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Adapted from LiteMol
 * Copyright (c) 2016 - now David Sehnal, licensed under Apache 2.0, See LICENSE file for more info.
 */

// import { throttle } from 'rxjs/operators';
// import { interval } from 'rxjs';

import { shallowClone } from 'mol-util';
import { Context } from '../../context/context'
import { Controller } from '../controller';

export const DefaultViewportOptions = {
    clearColor: { r: 1, g: 1, b: 1 },
    enableFog: true,
    cameraFOV: 30,
    cameraSpeed: 4
}
export type ViewportOptions = typeof DefaultViewportOptions

export class ViewportController extends Controller<ViewportOptions> {
    constructor(context: Context) {
        super(context, shallowClone(DefaultViewportOptions));
    }
}