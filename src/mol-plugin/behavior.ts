/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

export * from './behavior/behavior'
import * as Data from './behavior/built-in/data'
import * as Representation from './behavior/built-in/representation'

export const PluginBehaviors = {
    Data,
    Representation
}