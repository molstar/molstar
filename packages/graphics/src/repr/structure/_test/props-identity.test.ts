/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import fs from 'node:fs';
import path from 'node:path';
import { describe, expect, it } from 'vitest';
import { RuntimeContext } from '@molstar/core/task';
import { Color } from '@molstar/core/util/color';
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
import { SpacefillRepresentationProvider } from '../representation/spacefill.js';

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
  const parsed = await CIF.parse(fs.readFileSync(file, 'utf8').replace(/\r\n/g, '\n')).run();
  if (parsed.isError) throw new Error(parsed.message);
  const trajectory = await trajectoryFromMmCIF(parsed.result.blocks[0], parsed.result).run();
  return Structure.ofModel(trajectory.representative);
}

function themeOf(structure: Structure, color: number) {
  return Theme.create(
    registries,
    { structure },
    { colorTheme: { name: 'uniform', params: { value: Color(color) } }, sizeTheme: { name: 'uniform' } },
  );
}

describe('a representation update without prop changes', () => {
  it.each([
    ['units visuals (spacefill)', SpacefillRepresentationProvider, undefined],
    ['a complex visual (bonds per structure)', BallAndStickRepresentationProvider, ['structure-intra-bond']],
  ] as const)('keeps the values from the props of %s; a prop change after it applies', async (_, provider, visuals) => {
    const structure = await load1cbs();
    const defaults = PD.getDefaultValues(provider.getParams(registries, structure) as PD.Params);
    const params = visuals ? { ...defaults, visuals: [...visuals] } : defaults;
    const repr = (provider as typeof SpacefillRepresentationProvider).factory(
      { graphics: graphicsStub, ...registries },
      provider.getParams as typeof SpacefillRepresentationProvider.getParams,
    );
    repr.setTheme(themeOf(structure, 0xff0000));
    await repr.createOrUpdate(params as never, structure).runInContext(RuntimeContext.Synchronous);
    const objects = repr.renderObjects;
    expect(objects.length).toBeGreaterThan(0);
    // `uLod` is written by every `updateValues` (no compare): its version tells whether they ran
    const lodVersions = objects.map((o) => o.values.uLod.ref.version);

    // a recolor: the colors change, the values from the props stay
    repr.setTheme(themeOf(structure, 0x0000ff));
    await repr.createOrUpdate({}, structure).runInContext(RuntimeContext.Synchronous);
    expect(objects.map((o) => o.values.uLod.ref.version)).toEqual(lodVersions);
    for (const o of objects) expect(o.values.uColor.ref.value).toEqual(Color.toVec3Normalized([0, 0, 0], 0x0000ff));

    // then a prop change reaches every visual's values
    await repr.createOrUpdate({ alpha: 0.5 }, structure).runInContext(RuntimeContext.Synchronous);
    for (const o of repr.renderObjects) expect(o.values.uAlpha.ref.value).toBe(0.5);
  });
});
