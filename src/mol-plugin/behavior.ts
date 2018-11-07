/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

export * from './behavior/behavior'
import * as Data from './behavior/data'
import * as Representation from './behavior/representation'

export const PluginBehaviors = {
    Data,
    Representation
}