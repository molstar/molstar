/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { StateAction } from 'mol-state/action';
import { Transformer } from 'mol-state';
import { StateTransformParameters } from './ui/state/parameters';

export { PluginSpec }

interface PluginSpec {
    actions: PluginSpec.Action[],
    behaviors: PluginSpec.Behavior[]
}

namespace PluginSpec {
    export interface Action {
        action: StateAction | Transformer,
        customControl?: StateTransformParameters.Class,
        autoUpdate?: boolean
    }

    export function Action(action: StateAction | Transformer, params?: { customControl?: StateTransformParameters.Class, autoUpdate?: boolean }): Action {
        return { action, customControl: params && params.customControl, autoUpdate: params && params.autoUpdate };
    }

    export interface Behavior {
        transformer: Transformer,
        defaultParams?: any
    }

    export function Behavior<T extends Transformer>(transformer: T, defaultParams?: Transformer.Params<T>): Behavior {
        return { transformer, defaultParams };
    }
}