/**
 * Copyright (c) 2023-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Adam Midlik <midlik@gmail.com>
 * @author David Sehnal <david.sehnal@gmail.com>
 */

/**
 * Value palettes with Mol*'s named color lists: the evaluation (`value-palette-eval.ts`, which imports no lists),
 * re-exported, and the lists a palette's `colors` may name (`ValuePaletteLists`, MVS's `ColorListNameT`), which
 * `valuePalettePropsFromSpec` here reads. A consumer that names no Mol* list imports `value-palette-eval` and ships no
 * lists.
 */

import { Color } from './color.js';
import { ColorLists } from './lists.js';
import {
  type ValuePaletteProps,
  type ValuePalettePropsOptions,
  type ValuePaletteSpec,
  type ValuePaletteSpecDefaults,
  valuePalettePropsFromSpec as propsFromSpec,
} from './value-palette-eval.js';

export * from './value-palette-eval.js';

/** MVS's color list names (`ColorListNameT`) and the lists of `ColorLists` they name, in MVS's order. */
const ListNames = {
  // Sequential single-hue
  Reds: 'reds',
  Oranges: 'oranges',
  Greens: 'greens',
  Blues: 'blues',
  Purples: 'purples',
  Greys: 'greys',

  // Sequential multi-hue
  OrRd: 'orange-red',
  BuGn: 'blue-green',
  PuBuGn: 'purple-blue-green',
  GnBu: 'green-blue',
  PuBu: 'purple-blue',
  BuPu: 'blue-purple',
  RdPu: 'red-purple',
  PuRd: 'purple-red',
  YlOrRd: 'yellow-orange-red',
  YlOrBr: 'yellow-orange-brown',
  YlGn: 'yellow-green',
  YlGnBu: 'yellow-green-blue',

  Magma: 'magma',
  Inferno: 'inferno',
  Plasma: 'plasma',
  Viridis: 'viridis',
  Cividis: 'cividis',
  Turbo: 'turbo',
  Warm: 'warm',
  Cool: 'cool',
  CubehelixDefault: 'cubehelix-default',

  // Cyclical
  Rainbow: 'rainbow',
  Sinebow: 'sinebow',

  // Diverging
  RdBu: 'red-blue',
  RdGy: 'red-grey',
  PiYG: 'pink-yellow-green',
  BrBG: 'brown-white-green',
  PRGn: 'purple-green',
  PuOr: 'purple-orange',
  RdYlGn: 'red-yellow-green',
  RdYlBu: 'red-yellow-blue',
  Spectral: 'spectral',

  // Categorical
  Category10: 'category-10',
  Observable10: 'observable-10',
  Tableau10: 'tableau-10',

  Set1: 'set-1',
  Set2: 'set-2',
  Set3: 'set-3',
  Pastel1: 'pastel-1',
  Pastel2: 'pastel-2',
  Dark2: 'dark-2',
  Paired: 'paired',
  Accent: 'accent',

  // Additional lists, not standard for visualization in general, but commonly used for structures
  Chainbow: 'turbo-no-black',
} as const satisfies Record<string, keyof typeof ColorLists>;

export type ValuePaletteListName = keyof typeof ListNames;

const _lists: Partial<Record<ValuePaletteListName, readonly Color[]>> = {};

/**
 * The named color lists a palette's `colors` may name (MVS's `ColorListNameT`), as colors. Each list is converted when
 * first read (a getter per name), so loading the module converts none.
 */
export const ValuePaletteLists = {} as Readonly<Record<ValuePaletteListName, readonly Color[]>>;
for (const name of Object.keys(ListNames) as ValuePaletteListName[]) {
  Object.defineProperty(ValuePaletteLists, name, {
    enumerable: true,
    get: () => (_lists[name] ??= ColorLists[ListNames[name]].list.map(Color.fromColorListEntry)),
  });
}

/** A palette's props: `palette` over `defaults`; a `colors` name reads `ValuePaletteLists` unless `options.lists`. */
export function valuePalettePropsFromSpec(
  palette: ValuePaletteSpec | undefined | null,
  defaults: ValuePaletteSpecDefaults,
  options: ValuePalettePropsOptions = {},
): ValuePaletteProps {
  return propsFromSpec(palette, defaults, options.lists ? options : { ...options, lists: ValuePaletteLists });
}
