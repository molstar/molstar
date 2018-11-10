/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

export * from './behavior/behavior'

import * as State from './behavior/built-in/state'
import * as Representation from './behavior/built-in/representation'

export const BuiltInPluginBehaviors = {
    State,
}

export const PluginBehaviors = {
    Representation
}