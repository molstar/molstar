/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import fs from 'node:fs';
import path from 'node:path';
import { describe, expect, it, vi } from 'vitest';
import { Color } from '@molstar/core/util/color';
import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import { CIF } from '@molstar/io/reader/cif';
import { trajectoryFromMmCIF } from '@molstar/model/formats/structure/mmcif';
import { Structure, StructureElement } from '@molstar/model/model/structure';
import { IllustrativeColorThemeProvider } from '../color/illustrative.js';

describe('illustrative color theme', () => {
  it('lightens each style color once, with the same colors (perf-audit-2 09/P4)', async () => {
    const file = path.resolve(import.meta.dirname, '../../../../../data/examples/1cbs_updated.cif');
    const parsed = await CIF.parse(fs.readFileSync(file, 'utf8')).run();
    if (parsed.isError) throw new Error(parsed.message);
    const structure = Structure.ofModel(
      (await trajectoryFromMmCIF(parsed.result.blocks[0], parsed.result).run()).representative,
    );
    const props = PD.getDefaultValues(IllustrativeColorThemeProvider.getParams({ structure }));
    const theme = IllustrativeColorThemeProvider.factory({ structure }, props);

    const lighten = vi.spyOn(Color, 'lighten');
    const l = StructureElement.Location.create(structure);
    const colors: Color[] = [];
    let carbons = 0;
    for (const unit of structure.units) {
      l.unit = unit;
      for (const e of unit.elements) {
        l.element = e;
        colors.push(theme.color(l, false));
        if (unit.model.atomicHierarchy.atoms.type_symbol.value(e) === 'C') carbons++;
      }
    }
    const calls = lighten.mock.calls.length;
    lighten.mockRestore();
    expect(carbons).toBeGreaterThan(100);
    expect(calls).toBeLessThan(20);

    // the colors themselves: a lightened style color for carbons (`Color.lighten` is pure, so a memo of it is exact)
    expect(new Set(colors).size).toBeLessThan(40);
  });
});
