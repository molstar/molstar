/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

/** Index of an element in the model data */
export type ElementIndex = { readonly '@type': 'element-index' } & number
/** Index of a residue in the model data */
export type ResidueIndex = { readonly '@type': 'residue-index' } & number
/** Index of a chain in the model data */
export type ChainIndex = { readonly '@type': 'chain-index' } & number
/** Index of an entity in the model data */
export type EntityIndex = { readonly '@type': 'entity-index' } & number