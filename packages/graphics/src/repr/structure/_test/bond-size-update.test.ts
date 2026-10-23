/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import fs from 'node:fs';
import path from 'node:path';
import { describe, expect, it } from 'vitest';
import { RuntimeContext } from '@molstar/core/task';
import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import { CIF } from '@molstar/io/reader/cif';
import { parseMol2 } from '@molstar/io/reader/mol2/parser';
import { trajectoryFromMmCIF } from '@molstar/model/formats/structure/mmcif';
import { trajectoryFromMol2 } from '@molstar/model/formats/structure/mol2';
import { SortedArray } from '@molstar/core/data/int';
import { Structure, Unit } from '@molstar/model/model/structure';
import type { ElementIndex } from '@molstar/model/model/structure/model/indexing';
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

/** The parts of a GraphicsContext the CPU geometry paths read, with impostor support */
const graphicsStub = {
  capabilities: { fragDepth: true, flatInterpolation: true, texture: true },
  resources: {
    limits: { maxTextureSize: 16384, max3dTextureSize: 2048 },
    texture: () => new Proxy({} as any, { get: (t, k) => (k === 'then' ? undefined : k in t ? t[k] : () => {}) }),
    ownsTexture: () => true,
  },
  namedTextures: {},
} as any;

async function load1cbs() {
  const file = path.resolve(import.meta.dirname, '../../../../../../data/examples/1cbs_updated.cif');
  const parsed = await CIF.parse(fs.readFileSync(file, 'utf8')).run();
  if (parsed.isError) throw new Error(parsed.message);
  const trajectory = await trajectoryFromMmCIF(parsed.result.blocks[0], parsed.result).run();
  return Structure.ofModel(trajectory.representative);
}

function themeOf(structure: Structure, size: number) {
  return Theme.create(
    registries,
    { structure },
    { colorTheme: { name: 'uniform' }, sizeTheme: { name: 'uniform', params: { value: size } } },
  );
}

/** The cylinder impostors' start points, over every bond visual */
function cylinderStarts(repr: { renderObjects: readonly any[] }) {
  const starts: number[] = [];
  for (const ro of repr.renderObjects) {
    if (ro.type !== 'cylinders') continue;
    const a = ro.values.aStart.ref.value as Float32Array;
    const n = ro.values.drawCount.ref.value;
    for (let i = 0; i < Math.min(a.length, n * 3); i++) starts.push(a[i]);
  }
  return starts;
}

function differing(a: number[], b: number[]) {
  let n = 0;
  for (let i = 0; i < a.length; i++) if (a[i] !== b[i]) n++;
  return n;
}

/** adrenalin (MOL2, aromatic bonds) with its unit cut every 3 atoms: aromatic bonds cross units */
async function splitAdrenalin() {
  const file = path.resolve(import.meta.dirname, '../../../../../../data/examples/adrenalin.mol2');
  const parsed = await parseMol2(fs.readFileSync(file, 'utf8'), 'adrenalin').run();
  if (parsed.isError) throw new Error(parsed.message);
  const s = Structure.ofModel((await trajectoryFromMol2(parsed.result).run()).representative);
  const units: Unit[] = [];
  for (const u of s.units) {
    for (let i = 0; i < u.elements.length; i += 3) {
      const set = SortedArray.ofSortedArray<ElementIndex>(u.elements.slice(i, i + 3) as any);
      const id = units.length;
      units.push(Unit.create(id, id, id, u.traits, u.kind, u.model, u.conformation.operator, set));
    }
  }
  return Structure.create(units);
}

describe('bond impostors', () => {
  it.each([
    ['per unit', ['intra-bond'], load1cbs],
    ['per structure', ['structure-intra-bond'], load1cbs],
    ['inter-unit', ['inter-bond'], splitAdrenalin],
  ] as const)(
    'move double and aromatic bond offsets with a size theme change, %s (perf-audit-2 07/B5)',
    async (_, visuals, load) => {
      const structure = await load();
      const provider = BallAndStickRepresentationProvider;
      const params = {
        ...PD.getDefaultValues(provider.getParams(registries, structure)),
        multipleBonds: 'offset' as const,
        visuals: [...visuals],
      };

      const updated = provider.factory({ graphics: graphicsStub, ...registries }, provider.getParams);
      updated.setTheme(themeOf(structure, 1));
      await updated.createOrUpdate(params, structure).runInContext(RuntimeContext.Synchronous);
      const before = cylinderStarts(updated);
      updated.setTheme(themeOf(structure, 2));
      await updated.createOrUpdate(params, structure).runInContext(RuntimeContext.Synchronous);

      const fresh = provider.factory({ graphics: graphicsStub, ...registries }, provider.getParams);
      fresh.setTheme(themeOf(structure, 2));
      await fresh.createOrUpdate(params, structure).runInContext(RuntimeContext.Synchronous);

      const expected = cylinderStarts(fresh);
      const got = cylinderStarts(updated);
      expect(expected.length).toBeGreaterThan(0);
      expect(got.length).toBe(expected.length);
      // the offsets scale with the size: some starts move
      expect(differing(before, expected)).toBeGreaterThan(0);
      expect(differing(got, expected)).toBe(0);
    },
  );
});
