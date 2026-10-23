/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import fs from 'node:fs';
import path from 'node:path';
import { describe, expect, it } from 'vitest';
import { Mat4, Vec3 } from '@molstar/core/math/linear-algebra';
import { RuntimeContext } from '@molstar/core/task';
import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import { CIF } from '@molstar/io/reader/cif';
import { trajectoryFromMmCIF } from '@molstar/model/formats/structure/mmcif';
import { Structure } from '@molstar/model/model/structure';
import { ColorTheme } from '@molstar/graphics/theme/color';
import { BuiltInColorThemes } from '@molstar/graphics/theme/color/catalog';
import { SizeTheme } from '@molstar/graphics/theme/size';
import { BuiltInSizeThemes } from '@molstar/graphics/theme/size/catalog';
import { Theme } from '@molstar/graphics/theme/theme';
import { BallAndStickRepresentationProvider } from '../representation/ball-and-stick.js';

const registries = (() => {
  const colorThemeRegistry = ColorTheme.createRegistry();
  for (const p of Object.values(BuiltInColorThemes)) colorThemeRegistry.add(p as ColorTheme.Provider<any, any, any>);
  const sizeThemeRegistry = SizeTheme.createRegistry();
  for (const p of Object.values(BuiltInSizeThemes)) sizeThemeRegistry.add(p as SizeTheme.Provider<any, any>);
  return { colorThemeRegistry, sizeThemeRegistry };
})();

async function load1cbs() {
  const file = path.resolve(import.meta.dirname, '../../../../../../data/examples/1cbs_updated.cif');
  const parsed = await CIF.parse(fs.readFileSync(file, 'utf8')).run();
  if (parsed.isError) throw new Error(parsed.message);
  const trajectory = await trajectoryFromMmCIF(parsed.result.blocks[0], parsed.result).run();
  return Structure.ofModel(trajectory.representative);
}

describe('units representation transform', () => {
  it('applies the representation transform to visuals made after it was set (perf-audit-2 07/B2)', async () => {
    const structure = await load1cbs();
    expect(structure.units.length).toBeGreaterThan(1);
    const first = Structure.create([structure.units[0]], { parent: structure });

    const provider = BallAndStickRepresentationProvider;
    const repr = provider.factory({ ...registries }, provider.getParams);
    const themeProps = { colorTheme: { name: 'uniform' }, sizeTheme: { name: 'uniform' } };
    repr.setTheme(Theme.create(registries, { structure }, themeProps));
    const params = PD.getDefaultValues(provider.getParams(registries, structure));
    await repr.createOrUpdate(params, first).runInContext(RuntimeContext.Synchronous);

    const transform = Mat4.fromTranslation(Mat4(), Vec3.create(100, 0, 0));
    repr.setState({ transform });
    // a structure with more groups: visuals for the new groups are made now
    await repr.createOrUpdate(params, structure).runInContext(RuntimeContext.Synchronous);

    expect(repr.renderObjects.length).toBeGreaterThan(2);
    for (const ro of repr.renderObjects) {
      expect(Mat4.areEqual(ro.values.matrix.ref.value as Mat4, transform, 1e-6)).toBe(true);
    }
  });
});
