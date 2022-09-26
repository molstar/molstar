/**
 * Copyright (c) 2022 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

export interface wwPDBViewState {
    url: string,
    format: { name: string, isBinary: boolean },
    presetName?: string,
}