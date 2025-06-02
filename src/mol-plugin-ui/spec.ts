/**
 * Copyright (c) 2018-2025 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author Ventura Rivera <venturaxrivera@gmail.com>
 */


import { StateTransformParameters } from '../mol-plugin-ui/state/common';
import { CreateVolumeStreamingBehavior } from '../mol-plugin/behavior/dynamic/volume-streaming/transformers';
import { DefaultPluginSpec, PluginSpec } from '../mol-plugin/spec';
import { StateAction, StateTransformer } from '../mol-state';
import { VolumeStreamingCustomControls } from './custom/volume';
import { Loci } from '../mol-model/loci';
import { SequenceViewMode } from './sequence';

export interface PluginUISpec extends PluginSpec {
    customParamEditors?: [StateAction | StateTransformer, StateTransformParameters.Class][],
    components?: {
        controls?: PluginUISpec.LayoutControls
        remoteState?: 'none' | 'default',
        structureTools?: React.ComponentClass | React.FC,
        viewport?: {
            view?: React.ComponentClass | React.FC,
            controls?: React.ComponentClass | React.FC,
            snapshotDescription?: React.ComponentClass | React.FC,
        },
        sequenceViewer?: {
            view?: React.ComponentClass | React.FC
            modeOptions?: SequenceViewMode[],
            defaultMode?: SequenceViewMode,
        }
        hideTaskOverlay?: boolean,
        disableDragOverlay?: boolean,
        selectionTools?: {
            controls?: React.ComponentClass | React.FC,
            granularityOptions?: Loci.Granularity[],
            hide?: {
                granularity?: boolean,
                union?: boolean,
                subtract?: boolean,
                intersect?: boolean,
                set?: boolean,
                theme?: boolean,
                componentAdd?: boolean,
                componentRemove?: boolean,
                undo?: boolean,
                help?: boolean,
                cancel?: boolean,
            },
        },
    },
}

export namespace PluginUISpec {
    export interface LayoutControls {
        top?: React.ComponentClass | React.FC | 'none',
        left?: React.ComponentClass | React.FC | 'none',
        right?: React.ComponentClass | React.FC | 'none',
        bottom?: React.ComponentClass | React.FC | 'none'
    }
}

export const DefaultPluginUISpec = (): PluginUISpec => ({
    ...DefaultPluginSpec(),
    customParamEditors: [
        [CreateVolumeStreamingBehavior, VolumeStreamingCustomControls]
    ],
});