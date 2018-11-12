/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

export * from './behavior/behavior'

import * as StaticState from './behavior/static/state'
import * as StaticRepresentation from './behavior/static/representation'
import * as StaticCamera from './behavior/static/representation'

import * as DynamicRepresentation from './behavior/dynamic/representation'

export const BuiltInPluginBehaviors = {
    State: StaticState,
    Representation: StaticRepresentation,
    Camera: StaticCamera
}

export const PluginBehaviors = {
    Representation: DynamicRepresentation
}