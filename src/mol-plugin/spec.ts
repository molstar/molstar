/**
 * Copyright (c) 2018-2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { StateTransformer, StateAction } from '../mol-state';
import { StateTransformParameters } from '../mol-plugin-ui/state/common';
import { PluginLayoutStateProps } from './layout';
import { PluginStateAnimation } from '../mol-plugin-state/animation/model';
import { PluginConfigItem } from './config';

export { PluginSpec };

interface PluginSpec {
    actions: PluginSpec.Action[],
    behaviors: PluginSpec.Behavior[],
    animations?: PluginStateAnimation[],
    customParamEditors?: [StateAction | StateTransformer, StateTransformParameters.Class][],
    layout?: {
        initial?: Partial<PluginLayoutStateProps>,
        controls?: PluginSpec.LayoutControls
    },
    components?: {
        remoteState?: 'none' | 'default',
        structureTools?: React.ComponentClass,
        viewport?: {
            view?: React.ComponentClass,
            controls?: React.ComponentClass
        }
    },
    config?: [PluginConfigItem, unknown][]
}

namespace PluginSpec {
    export interface Action {
        action: StateAction | StateTransformer,
        customControl?: StateTransformParameters.Class,
        autoUpdate?: boolean
    }

    export function Action(action: StateAction | StateTransformer, params?: { customControl?: StateTransformParameters.Class, autoUpdate?: boolean }): Action {
        return { action, customControl: params && params.customControl, autoUpdate: params && params.autoUpdate };
    }

    export interface Behavior {
        transformer: StateTransformer,
        defaultParams?: any
    }

    export function Behavior<T extends StateTransformer>(transformer: T, defaultParams: Partial<StateTransformer.Params<T>> = {}): Behavior {
        return { transformer, defaultParams };
    }

    export interface LayoutControls {
        top?: React.ComponentClass | 'none',
        left?: React.ComponentClass | 'none',
        right?: React.ComponentClass | 'none',
        bottom?: React.ComponentClass | 'none'
    }
}