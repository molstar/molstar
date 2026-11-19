/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { readFileSync } from 'node:fs';
import { dirname, join, resolve } from 'node:path';
import { fileURLToPath } from 'node:url';
import { describe, expect, it } from 'vitest';
import { Color } from '../color.js';
import { ValuePaletteLists, valuePalettePropsFromSpec } from '../value-palette.js';
import {
  type ValuePaletteSpecDefaults,
  valuePalettePropsFromSpec as evalPropsFromSpec,
} from '../value-palette-eval.js';

const here = dirname(fileURLToPath(import.meta.url));

/** The relative modules a source file imports, transitively (static `import`/`export ... from` only). */
function closure(file: string, seen = new Set<string>()): Set<string> {
  if (seen.has(file)) return seen;
  seen.add(file);
  const src = readFileSync(file, 'utf8');
  for (const m of src.matchAll(/(?:import|export)\s[^;]*?from\s+'(\.[^']+)'/g)) {
    closure(resolve(dirname(file), m[1].replace(/\.js$/, '.ts')), seen);
  }
  return seen;
}

const Defaults: ValuePaletteSpecDefaults = {
  categorical: {
    kind: 'categorical',
    colors: [],
    repeat_color_list: false,
    sort: 'none',
    sort_direction: 'ascending',
    case_insensitive: false,
    missing_color: null,
  },
  discrete: { kind: 'discrete', colors: [], reverse: false, mode: 'normalized', value_domain: [null, null] },
  continuous: {
    kind: 'continuous',
    colors: [],
    reverse: false,
    mode: 'normalized',
    value_domain: [0, 1],
    underflow_color: null,
    overflow_color: null,
  },
};

describe('value-palette-eval', () => {
  it('imports no color lists (a consumer of the evaluation ships none)', () => {
    const files = [...closure(join(here, '../value-palette-eval.ts'))].map((f) => f.replace(/\\/g, '/'));
    expect(files.some((f) => f.endsWith('/util/color/value-palette-eval.ts'))).toBe(true);
    expect(files.filter((f) => f.endsWith('/util/color/lists.ts'))).toEqual([]);
  });

  it('reads a list name only from the lists it is given', () => {
    const unknown: string[] = [];
    const props = evalPropsFromSpec({ kind: 'continuous', colors: 'Viridis' }, Defaults, {
      onUnknownName: (n) => unknown.push(n),
    });
    expect(unknown).toEqual(['Viridis']);
    expect(props.name === 'continuous' && props.params.colors).toEqual([]);
    const given = evalPropsFromSpec({ kind: 'continuous', colors: 'Viridis' }, Defaults, { lists: ValuePaletteLists });
    expect(given).toEqual(valuePalettePropsFromSpec({ kind: 'continuous', colors: 'Viridis' }, Defaults));
  });

  it('value-palette names the lists by default, and a given `lists` wins', () => {
    const p = valuePalettePropsFromSpec({ kind: 'categorical', colors: 'Set1' }, Defaults);
    expect(p.name === 'categorical' && p.params.colors.params.length).toBe(ValuePaletteLists.Set1.length);
    const mine = valuePalettePropsFromSpec({ kind: 'categorical', colors: 'Mine' }, Defaults, {
      lists: { Mine: [Color(0xff0000)] },
    });
    expect(mine.name === 'categorical' && mine.params.colors.params.length).toBe(1);
  });
});
