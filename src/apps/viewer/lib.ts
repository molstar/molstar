/**
 * Copyright (c) 2025 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import * as Structure from '../../mol-model/structure';
import { DataLoci, EveryLoci, Loci } from '../../mol-model/loci';
import { Volume } from '../../mol-model/volume';
import { Shape, ShapeGroup } from '../../mol-model/shape';
import * as LinearAlgebra3D from '../../mol-math/linear-algebra/3d';
import { PluginContext } from '../../mol-plugin/context';
import { PluginConfig } from '../../mol-plugin/config';
import { PluginBehavior } from '../../mol-plugin/behavior';
import { DefaultPluginSpec, PluginSpec } from '../../mol-plugin/spec';
import { DefaultPluginUISpec } from '../../mol-plugin-ui/spec';
import { PluginStateObject, PluginStateTransform } from '../../mol-plugin-state/objects';
import { StateTransforms } from '../../mol-plugin-state/transforms';
import { StateActions } from '../../mol-plugin-state/actions';
import { PluginExtensions } from './extensions';

export const lib = {
    structure: {
        ...Structure,
    },
    volume: {
        Volume,
    },
    shape: {
        Shape,
        ShapeGroup,
    },
    loci: {
        Loci,
        DataLoci,
        EveryLoci,
    },
    math: {
        LinearAlgebra: {
            ...LinearAlgebra3D,
        }
    },
    plugin: {
        PluginContext,
        PluginConfig,
        PluginBehavior,
        PluginSpec,
        PluginStateObject,
        PluginStateTransform,
        StateTransforms,
        StateActions,
        DefaultPluginSpec,
        DefaultPluginUISpec,
    },
    extensions: {
        ...PluginExtensions
    }
};