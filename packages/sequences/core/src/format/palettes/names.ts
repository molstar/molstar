/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

/**
 * Mol*'s color names (X11 names, http://www.w3.org/TR/css3-color/#svg-color, with Mol*'s additions), the names a
 * document `Color` may use: `@molstar/core/util/color/names`.
 */

import { ColorNames as MolstarColorNames } from '@molstar/core/util/color/names';

/** Lower-case name → 0xRRGGBB. */
export const ColorNames: ReadonlyMap<string, number> = new Map(Object.entries(MolstarColorNames));
