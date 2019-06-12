/**
 * Copyright (c) 2018-2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { ConsoleLogger } from '../../mol-util/console-logger';

export function getParam<T>(params: any, ...path: string[]): T | undefined {
    try {
        let current = params;
        for (const p of path) {
            if (typeof current === 'undefined') return;
            current = current[p];
        }
        return current;
    } catch (e) {
        ConsoleLogger.error('Config', `Unable to retrieve property ${path.join('.')} from ${JSON.stringify(params)}`);
    }
}