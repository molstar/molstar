/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import fs from 'node:fs';
import path from 'node:path';
import { describe, expect, it, vi } from 'vitest';
import { hashFnv32a } from '@molstar/core/data/util/hash-functions';
import type { PositionData } from '@molstar/core/math/geometry';
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
import { MolecularSurfaceMeshParams } from '../visual/molecular-surface-mesh.js';
import { computeUnitMolecularSurface } from '../visual/util/molecular-surface.js';

const positions: PositionData[] = [];
vi.mock('@molstar/core/math/geometry/molecular-surface', async (importOriginal) => {
  const original = await importOriginal<typeof import('@molstar/core/math/geometry/molecular-surface')>();
  return {
    ...original,
    calcMolecularSurface: (ctx: any, position: PositionData, ...rest: any[]) => {
      positions.push(position);
      return (original.calcMolecularSurface as any)(ctx, position, ...rest);
    },
  };
});

async function load1cbs() {
  const file = path.resolve(import.meta.dirname, '../../../../../../data/examples/1cbs_updated.cif');
  const parsed = await CIF.parse(fs.readFileSync(file, 'utf8')).run();
  if (parsed.isError) throw new Error(parsed.message);
  const trajectory = await trajectoryFromMmCIF(parsed.result.blocks[0], parsed.result).run();
  return Structure.ofModel(trajectory.representative);
}

function hash(a: ArrayLike<number>) {
  const f = Float32Array.from(a);
  return (hashFnv32a(new Uint8Array(f.buffer)) >>> 0).toString(16);
}

describe('unit molecular surface', () => {
  it('gives the surface per-unit arrays, not ones sized to the model (perf-audit-2 07/P6)', async () => {
    const structure = await load1cbs();
    const colorThemeRegistry = ColorTheme.createRegistry();
    for (const p of Object.values(BuiltInColorThemes)) colorThemeRegistry.add(p as ColorTheme.Provider<any, any, any>);
    const sizeThemeRegistry = SizeTheme.createRegistry();
    for (const p of Object.values(BuiltInSizeThemes)) sizeThemeRegistry.add(p as SizeTheme.Provider<any, any>);
    const theme = Theme.create(
      { colorThemeRegistry, sizeThemeRegistry },
      { structure },
      { colorTheme: { name: 'uniform' }, sizeTheme: { name: 'physical' } },
    );
    const sizeTheme = theme.size;
    const props = { ...PD.getDefaultValues(MolecularSurfaceMeshParams), resolution: 1 };
    const fields: string[] = [];
    for (const unit of structure.units.slice(0, 2)) {
      positions.length = 0;
      const density = await computeUnitMolecularSurface(structure, unit, sizeTheme, props).runInContext(
        RuntimeContext.Synchronous,
      );
      expect(positions).toHaveLength(1);
      expect(positions[0].radius!.length).toBe(unit.elements.length);
      expect(positions[0].x.length).toBe(unit.elements.length);
      fields.push(`${hash(density.field.space.dimensions)} ${hash(density.field.data)} ${hash(density.idField.data)}`);
    }
    // the fields as before (pinned on the model-sized arrays: the same hashes)
    expect(fields).toEqual(['ffcb98b3 377ad67b ea0c6c59', '738ddfb2 3f3bf7ea 65de0692']);
  });
});
