/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import fs from 'node:fs';
import path from 'node:path';
import { describe, expect, it, vi } from 'vitest';
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
import { CartoonRepresentationProvider } from '../representation/cartoon.js';
import { SpacefillRepresentationProvider } from '../representation/spacefill.js';

// counts the walks of the visuals' foreign-resource checks
const checks = vi.hoisted(() => ({ count: 0 }));
vi.mock('@molstar/graphics/gl/render-object', async (importOriginal) => {
  const original = await importOriginal<typeof import('@molstar/graphics/gl/render-object')>();
  return {
    ...original,
    createForeignResourcesCheck: () => {
      const check = original.createForeignResourcesCheck();
      return (o: any, graphics: any) => (checks.count++, check(o, graphics));
    },
  };
});

const registries = (() => {
  const colorThemeRegistry = ColorTheme.createRegistry();
  for (const p of Object.values(BuiltInColorThemes)) colorThemeRegistry.add(p as ColorTheme.Provider<any, any, any>);
  const sizeThemeRegistry = SizeTheme.createRegistry();
  for (const p of Object.values(BuiltInSizeThemes)) sizeThemeRegistry.add(p as SizeTheme.Provider<any, any>);
  return { colorThemeRegistry, sizeThemeRegistry };
})();

/** The parts of a GraphicsContext the CPU geometry paths read */
const graphicsStub = {
  capabilities: { fragDepth: true, flatInterpolation: true, texture: true },
  resources: {
    limits: { maxTextureSize: 16384, max3dTextureSize: 2048 },
    texture: () => new Proxy({} as any, { get: (t, k) => (k === 'then' ? undefined : k in t ? t[k] : () => {}) }),
    ownsTexture: () => true,
  },
  namedTextures: {},
} as any;

describe('foreign resources', () => {
  it.each([
    ['units', CartoonRepresentationProvider],
    ['complex', SpacefillRepresentationProvider],
  ] as const)('the first update after a build checks no %s visual (perf-audit-2 07/13/P4)', async (_, provider) => {
    const file = path.resolve(import.meta.dirname, '../../../../../../data/examples/1cbs_updated.cif');
    const parsed = await CIF.parse(fs.readFileSync(file, 'utf8')).run();
    if (parsed.isError) throw new Error(parsed.message);
    const structure = Structure.ofModel(
      (await trajectoryFromMmCIF(parsed.result.blocks[0], parsed.result).run()).representative,
    );
    const params: any = PD.getDefaultValues(provider.getParams(registries, structure));
    if (provider === SpacefillRepresentationProvider) params.visuals = ['structure-element-sphere'];
    const theme = (name: string) =>
      Theme.create(registries, { structure }, { colorTheme: { name }, sizeTheme: { name: 'uniform' } });

    const repr = provider.factory({ graphics: graphicsStub, ...registries }, provider.getParams);
    repr.setTheme(theme('uniform'));
    await repr.createOrUpdate(params, structure).runInContext(RuntimeContext.Synchronous);
    expect(repr.renderObjects.length).toBeGreaterThan(0);

    checks.count = 0;
    repr.setTheme(theme('chain-id'));
    await repr.createOrUpdate(params, structure).runInContext(RuntimeContext.Synchronous);
    expect(checks.count).toBe(0);

    // another context is checked
    const other = { ...graphicsStub };
    repr.setGraphics?.(other);
    await repr.createOrUpdate(params, structure).runInContext(RuntimeContext.Synchronous);
    expect(checks.count).toBeGreaterThan(0);
  });
});
