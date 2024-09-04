/**
 * Copyright (c) 2024 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Adam Midlik <midlik@gmail.com>
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Aliaksei Chareshneu <chareshneu.tech@gmail.com>
 */

import { PluginCommands } from '../../../mol-plugin/commands';
import { PluginContext } from '../../../mol-plugin/context';


class _MVSErrorContext {
    errors: string[];

    add(error: string) {
        this.errors.push(error);
    }
    reset() {
        // Reset errors
        this.errors = [];
    }
    report(ctx: PluginContext) {
        for (const error of this.errors) {
            ctx.log.warn(error);
            PluginCommands.Toast.Show(ctx, {
                title: 'Error',
                message: error,
                timeoutMs: 10000
            });
        }

    }
}

export const MVSErrorContext = new _MVSErrorContext();
