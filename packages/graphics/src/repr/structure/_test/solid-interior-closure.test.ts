/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Gianluca Tomasello <giagitom@gmail.com>
 */

import fs from 'node:fs';
import path from 'node:path';
import { describe, expect, it } from 'vitest';
import { RuntimeContext } from '@molstar/core/task';
import { DefaultAssetManager } from '@molstar/core/util/assets';
import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import { CIF } from '@molstar/io/reader/cif';
import { trajectoryFromMmCIF } from '@molstar/model/formats/structure/mmcif';
import { Structure, StructureElement, StructureProperties } from '@molstar/model/model/structure';
import type { CustomProperty } from '@molstar/model/props/common/custom-property';
import { ColorTheme } from '@molstar/graphics/theme/color';
import { BuiltInColorThemes } from '@molstar/graphics/theme/color/catalog';
import { SizeTheme } from '@molstar/graphics/theme/size';
import { BuiltInSizeThemes } from '@molstar/graphics/theme/size/catalog';
import { Theme } from '@molstar/graphics/theme/theme';
import type { MeshValues } from '@molstar/graphics/gl/renderable/mesh';
import type { Representation } from '../../representation.js';
import { BlobSurfaceRepresentationProvider } from '../representation/blob-surface.js';

const registries = (() => {
  const colorThemeRegistry = ColorTheme.createRegistry();
  for (const p of Object.values(BuiltInColorThemes)) colorThemeRegistry.add(p as ColorTheme.Provider<any, any, any>);
  const sizeThemeRegistry = SizeTheme.createRegistry();
  for (const p of Object.values(BuiltInSizeThemes)) sizeThemeRegistry.add(p as SizeTheme.Provider<any, any>);
  return { colorThemeRegistry, sizeThemeRegistry };
})();
const propCtx: CustomProperty.Context = {
  runtime: RuntimeContext.Synchronous,
  assetManager: new DefaultAssetManager(),
};

async function load1cbs() {
  const file = path.resolve(import.meta.dirname, '../../../../../../data/examples/1cbs_updated.cif');
  const parsed = await CIF.parse(fs.readFileSync(file, 'utf8')).run();
  if (parsed.isError) throw new Error(parsed.message);
  const trajectory = await trajectoryFromMmCIF(parsed.result.blocks[0], parsed.result).run();
  return Structure.ofModel(trajectory.representative);
}

/** the residues from `min` to `max` of `parent`, with `parent` as its parent (as the plugin's components) */
function residues(parent: Structure, min: number, max: number) {
  const seqId = StructureProperties.residue.label_seq_id;
  const l = StructureElement.Location.create(parent);
  const elements: { unit: (typeof parent.units)[number]; indices: StructureElement.UnitIndex[] }[] = [];
  for (const unit of parent.units) {
    l.unit = unit;
    const indices: StructureElement.UnitIndex[] = [];
    for (let i = 0; i < unit.elements.length; ++i) {
      l.element = unit.elements[i];
      const s = seqId(l);
      if (s >= min && s <= max) indices.push(i as StructureElement.UnitIndex);
    }
    if (indices.length) elements.push({ unit, indices });
  }
  const loci = StructureElement.Loci(parent, elements as any);
  return Structure.create(StructureElement.Loci.toStructure(loci).units, { parent });
}

/** the meshes of the representation made with `props`, then updated with `update` when given */
async function meshesOf(provider: any, structure: Structure, props: object = {}, update?: object) {
  await provider.ensureCustomProperties?.attach(propCtx, structure);
  const themeProps = {
    colorTheme: { name: 'uniform', params: {} },
    sizeTheme: { name: provider.defaultSizeTheme.name, params: provider.defaultSizeTheme.props ?? {} },
  };
  await Theme.ensureDependencies(propCtx, registries, { structure }, themeProps);
  const repr = provider.factory({ ...registries }, provider.getParams) as Representation<Structure>;
  repr.setTheme(Theme.create(registries, { structure }, themeProps));
  const params = PD.getDefaultValues(provider.getParams(registries, structure));
  await repr.createOrUpdate({ ...params, ...props }, structure).run();
  if (update) await repr.createOrUpdate({ ...params, ...props, ...update }, structure).run();
  const meshes = repr.renderObjects.filter((o) => o.type === 'mesh').map((o) => o.values as MeshValues);
  repr.destroy();
  return meshes;
}

/** edges by vertex position (zero-length triangles left out): those not shared by two triangles, and directed edges
 * without their reverse */
function openEdges(v: MeshValues) {
  const vb = v.aPosition.ref.value;
  const ib = v.elements.ref.value;
  const weld = new Int32Array(v.uVertexCount.ref.value);
  const keys = new Map<string, number>();
  for (let i = 0; i < weld.length; ++i) {
    const k = `${vb[i * 3].toFixed(3)},${vb[i * 3 + 1].toFixed(3)},${vb[i * 3 + 2].toFixed(3)}`;
    if (!keys.has(k)) keys.set(k, keys.size);
    weld[i] = keys.get(k)!;
  }
  const counts = new Map<number, number>();
  for (let t = 0, tl = v.drawCount.ref.value; t < tl; t += 3) {
    const tri = [weld[ib[t]], weld[ib[t + 1]], weld[ib[t + 2]]];
    if (tri[0] === tri[1] || tri[1] === tri[2] || tri[2] === tri[0]) continue;
    for (let j = 0; j < 3; ++j) {
      const k = tri[j] * 2 ** 26 + tri[(j + 1) % 3];
      counts.set(k, (counts.get(k) ?? 0) + 1);
    }
  }
  let unpaired = 0;
  for (const [k, c] of counts) {
    const a = Math.floor(k / 2 ** 26);
    const b = k % 2 ** 26;
    if ((counts.get(b * 2 ** 26 + a) ?? 0) !== c) ++unpaired;
  }
  return unpaired;
}

describe('blob surface with includeParent', async () => {
  const root = await load1cbs();
  const cut = residues(root, 20, 60);

  it('is closed against the parent', async () => {
    const [mesh] = await meshesOf(BlobSurfaceRepresentationProvider, cut, { includeParent: true, blobSize: 8 });
    expect(mesh.uClosureVertexCount.ref.value).toBeGreaterThan(0);
    expect(openEdges(mesh)).toBe(0);
  });
});
