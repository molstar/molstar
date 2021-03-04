/**
 * Copyright (c) 2018-2021 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { PartialCanvas3DProps } from '../mol-canvas3d/canvas3d';
import { AnimateAssemblyUnwind } from '../mol-plugin-state/animation/built-in/assembly-unwind';
import { AnimateCameraSpin } from '../mol-plugin-state/animation/built-in/camera-spin';
import { AnimateModelIndex } from '../mol-plugin-state/animation/built-in/model-index';
import { AnimateStateSnapshots } from '../mol-plugin-state/animation/built-in/state-snapshots';
import { PluginStateAnimation } from '../mol-plugin-state/animation/model';
import { DataFormatProvider } from '../mol-plugin-state/formats/provider';
import { StateTransformer } from '../mol-state';
import { PluginBehaviors } from './behavior';
import { StructureFocusRepresentation } from './behavior/dynamic/selection/structure-focus-representation';
import { PluginConfigItem } from './config';
import { PluginLayoutStateProps } from './layout';

export { PluginSpec };

interface PluginSpec {
    behaviors: PluginSpec.Behavior[],
    animations?: PluginStateAnimation[],
    customFormats?: [string, DataFormatProvider][],
    canvas3d?: PartialCanvas3DProps,
    layout?: {
        initial?: Partial<PluginLayoutStateProps>,
    },
    config?: [PluginConfigItem, unknown][]
}

namespace PluginSpec {
    export interface Behavior {
        transformer: StateTransformer,
        defaultParams?: any
    }

    export function Behavior<T extends StateTransformer>(transformer: T, defaultParams: Partial<StateTransformer.Params<T>> = {}): Behavior {
        return { transformer, defaultParams };
    }
}

export const DefaultPluginSpec = (): PluginSpec => ({
    behaviors: [
        PluginSpec.Behavior(PluginBehaviors.Representation.HighlightLoci),
        PluginSpec.Behavior(PluginBehaviors.Representation.SelectLoci),
        PluginSpec.Behavior(PluginBehaviors.Representation.DefaultLociLabelProvider),
        PluginSpec.Behavior(PluginBehaviors.Representation.FocusLoci),
        PluginSpec.Behavior(PluginBehaviors.Camera.FocusLoci),
        PluginSpec.Behavior(StructureFocusRepresentation),

        PluginSpec.Behavior(PluginBehaviors.CustomProps.StructureInfo),
        PluginSpec.Behavior(PluginBehaviors.CustomProps.AccessibleSurfaceArea),
        PluginSpec.Behavior(PluginBehaviors.CustomProps.Interactions),
        PluginSpec.Behavior(PluginBehaviors.CustomProps.SecondaryStructure),
        PluginSpec.Behavior(PluginBehaviors.CustomProps.ValenceModel),
        PluginSpec.Behavior(PluginBehaviors.CustomProps.CrossLinkRestraint),
    ],
    animations: [
        AnimateModelIndex,
        AnimateCameraSpin,
        AnimateStateSnapshots,
        AnimateAssemblyUnwind
    ]
});