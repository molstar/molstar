/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

/**
 * MVS's named color lists (`ColorListNameT`), keyed by MVS's names: Mol*'s (`@molstar/core/util/color/value-palette`,
 * `ValuePaletteLists`), which MVS's runtime reads too.
 */

import { ValuePaletteLists } from '@molstar/core/util/color/value-palette';
import type { Color } from './color.js';

export const ColorLists: Readonly<Record<keyof typeof ValuePaletteLists, readonly Color[]>> = ValuePaletteLists;

export type ColorListName = keyof typeof ColorLists;
export const ColorListNames = Object.keys(ColorLists) as readonly ColorListName[];
