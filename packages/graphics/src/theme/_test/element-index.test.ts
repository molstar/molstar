/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import fs from 'node:fs';
import path from 'node:path';
import { describe, expect, it } from 'vitest';
import { OrderedSet, SortedArray } from '@molstar/core/data/int';
import { getPalette } from '@molstar/core/util/color/palette';
import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import { CIF } from '@molstar/io/reader/cif';
import { trajectoryFromMmCIF } from '@molstar/model/formats/structure/mmcif';
import { Structure, StructureElement, Unit } from '@molstar/model/model/structure';
import type { ElementIndex } from '@molstar/model/model/structure/model/indexing';
import { ElementIndexColorTheme, ElementIndexColorThemeParams } from '../color/element-index.js';

async function load1cbs() {
  const file = path.resolve(import.meta.dirname, '../../../../../data/examples/1cbs_updated.cif');
  const parsed = await CIF.parse(fs.readFileSync(file, 'utf8').replace(/\r\n/g, '\n')).run();
  if (parsed.isError) throw new Error(parsed.message);
  const trajectory = await trajectoryFromMmCIF(parsed.result.blocks[0], parsed.result).run();
  return Structure.ofModel(trajectory.representative);
}

/** `s` with its units' every other element: elements with gaps */
function withGaps(s: Structure) {
  const units = s.units.map((u, i) => {
    const elements = SortedArray.ofSortedArray<ElementIndex>(Array.from(u.elements).filter((_, j) => j % 2 === 0));
    return Unit.create(i, i, i, u.traits, u.kind, u.model, u.conformation.operator, elements);
  });
  return Structure.create(units);
}

describe('element index colors', () => {
  it('color each element by its index in the structure, units of consecutive elements and of gaps alike', async () => {
    const full = await load1cbs();
    for (const structure of [full, withGaps(full)]) {
      const props = PD.getDefaultValues(ElementIndexColorThemeParams);
      const theme = ElementIndexColorTheme({ structure }, props);
      const { units } = structure.root;
      let elementCount = 0;
      for (const u of units) elementCount += u.elements.length;
      const palette = getPalette(elementCount, props);

      const location = StructureElement.Location.create(structure);
      let base = 0;
      let compared = 0;
      for (const unit of units) {
        location.unit = unit;
        for (let j = 0; j < unit.elements.length; ++j) {
          location.element = unit.elements[j];
          const expected = palette.color(base + OrderedSet.findPredecessorIndex(unit.elements, location.element));
          if (theme.color(location, false) !== expected) throw new Error(`unit ${unit.id}, element ${j}`);
          compared += 1;
        }
        base += unit.elements.length;
      }
      expect(compared).toBe(elementCount);
    }
  });
});
