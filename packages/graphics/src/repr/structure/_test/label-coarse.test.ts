/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import fs from 'node:fs';
import path from 'node:path';
import { afterAll, beforeAll, describe, expect, it } from 'vitest';
import { RuntimeContext } from '@molstar/core/task';
import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import { CIF } from '@molstar/io/reader/cif';
import { trajectoryFromMmCIF } from '@molstar/model/formats/structure/mmcif';
import { Structure, Unit } from '@molstar/model/model/structure';
import { clearFontAtlasCache, setCanvasModule } from '@molstar/graphics/geo/geometry/text/font-atlas';
import { ColorTheme } from '@molstar/graphics/theme/color';
import { BuiltInColorThemes } from '@molstar/graphics/theme/color/catalog';
import { SizeTheme } from '@molstar/graphics/theme/size';
import { BuiltInSizeThemes } from '@molstar/graphics/theme/size/catalog';
import { Theme } from '@molstar/graphics/theme/theme';
import { LabelRepresentationProvider } from '../representation/label.js';

// a font that draws nothing: the labels are counted, not looked at
beforeAll(() => {
  const context = {
    measureText: (text: string) => ({ width: 10 * text.length, actualBoundingBoxLeft: 0, actualBoundingBoxRight: 10 }),
    clearRect() {},
    fillText() {},
    getImageData: (_x: number, _y: number, w: number, h: number) => ({ data: new Uint8ClampedArray(w * h * 4) }),
  };
  setCanvasModule({ createCanvas: () => ({ getContext: () => context }) });
});
afterAll(() => {
  setCanvasModule(undefined);
  clearFontAtlasCache();
});

const registries = (() => {
  const colorThemeRegistry = ColorTheme.createRegistry();
  for (const p of Object.values(BuiltInColorThemes)) colorThemeRegistry.add(p as ColorTheme.Provider<any, any, any>);
  const sizeThemeRegistry = SizeTheme.createRegistry();
  for (const p of Object.values(BuiltInSizeThemes)) sizeThemeRegistry.add(p as SizeTheme.Provider<any, any>);
  return { colorThemeRegistry, sizeThemeRegistry };
})();

async function load8zzh() {
  const file = path.resolve(import.meta.dirname, '../../../../../../data/examples/ihm-coarse-8zzh.cif');
  const parsed = await CIF.parse(fs.readFileSync(file, 'utf8')).run();
  if (parsed.isError) throw new Error(parsed.message);
  const trajectory = await trajectoryFromMmCIF(parsed.result.blocks[0], parsed.result).run();
  return Structure.ofModel(trajectory.representative);
}

/** The characters drawn by a label representation of `structure` at `level`. */
async function labelChars(structure: Structure, level: 'chain' | 'residue' | 'element') {
  const provider = LabelRepresentationProvider;
  const repr = provider.factory({ ...registries }, provider.getParams);
  repr.setTheme(
    Theme.create(registries, { structure }, { colorTheme: { name: 'uniform' }, sizeTheme: { name: 'uniform' } }),
  );
  const params = { ...PD.getDefaultValues(provider.getParams(registries, structure)), level };
  await repr.createOrUpdate(params, structure).runInContext(RuntimeContext.Synchronous);
  let chars = 0;
  for (const ro of repr.renderObjects) chars += (ro.values.drawCount.ref.value as number) / 6;
  return chars;
}

describe('label representation on coarse units', () => {
  it('labels coarse chains by their asym id and skips them at residue and element level', async () => {
    const structure = await load8zzh();
    expect(structure.units.every((u) => !Unit.isAtomic(u))).toBe(true);
    expect(await labelChars(structure, 'chain')).toBeGreaterThan(0);
    expect(await labelChars(structure, 'residue')).toBe(0);
    expect(await labelChars(structure, 'element')).toBe(0);
  });
});
