/**
 * Copyright (c) 2024 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Adam Midlik <midlik@gmail.com>
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Aliaksei Chareshneu <chareshneu.tech@gmail.com>
 */

import { PluginContext } from '../../../mol-plugin/context';


class _MVSErrorContext {
    add(error: any) {

    }
    reset() {
        // Reset errors
    }
    report(ctx: PluginContext) {

    }
}

export const MVSErrorContext = new _MVSErrorContext();
