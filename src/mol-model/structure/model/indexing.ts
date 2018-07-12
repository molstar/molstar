/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

export type ElementIndex = { readonly '@type': 'element-index' } & number
export type ResidueIndex = { readonly '@type': 'residue-index' } & number
export type ChainIndex = { readonly '@type': 'chain-index' } & number
export type EntityIndex = { readonly '@type': 'entity-index' } & number