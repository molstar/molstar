/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import fs from 'node:fs';
import path from 'node:path';
import { describe, expect, it } from 'vitest';
import { hashFnv32a } from '@molstar/core/data/util/hash-functions';
import { Vec3 } from '@molstar/core/math/linear-algebra';
import { SyncRuntimeContext } from '@molstar/core/task/execution/synchronous';
import { AssetManager } from '@molstar/core/util/assets';
import { MarkerAction } from '@molstar/core/util/marker-action';
import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import { CIF } from '@molstar/io/reader/cif';
import { trajectoryFromMmCIF } from '@molstar/model/formats/structure/mmcif';
import type { Loci } from '@molstar/model/model/loci';
import { Structure, StructureElement } from '@molstar/model/model/structure';
import type { Model } from '@molstar/model/model/structure/model';
import { StructureSymmetry } from '@molstar/model/model/structure/structure/symmetry';
import { InteractionsProvider } from '@molstar/model/props/computed/interactions';
import { Bridges, Interactions } from '@molstar/model/props/computed/interactions/interactions';
import { BridgeProviderParams, ContactProviderParams } from '@molstar/model/props/computed/interactions/params';
import type { GraphicsRenderObject } from '@molstar/graphics/gl/render-object';
import type { Representation } from '@molstar/graphics/repr/representation';
import { ColorTheme } from '@molstar/graphics/theme/color';
import { BuiltInColorThemes } from '@molstar/graphics/theme/color/catalog';
import { SizeTheme } from '@molstar/graphics/theme/size';
import { BuiltInSizeThemes } from '@molstar/graphics/theme/size/catalog';
import { Theme } from '@molstar/graphics/theme/theme';
import { InteractionTypeColorThemeProvider } from '../../themes/interaction-type.js';
import { InteractionsRepresentationProvider } from '../interactions.js';

// Pins what the interactions representation draws and how its groups map to loci (picking) and back (marking), so
// that the cylinder impostor visuals and the mesh fallback can be checked against each other (09-graphics-theme-props
// P1).

const examples = path.resolve(import.meta.dirname, '../../../../../../../data/examples');

async function loadCif(name: string): Promise<Model> {
  const parsed = await CIF.parse(fs.readFileSync(path.join(examples, name), 'utf8')).run();
  if (parsed.isError) throw new Error(parsed.message);
  return (await trajectoryFromMmCIF(parsed.result.blocks[0], parsed.result).run()).representative;
}

/** Every contact and bridge provider on, with its default params. */
function allOn(params: Record<string, PD.Mapped<any>>) {
  const ret: Record<string, { name: 'on'; params: any }> = {};
  for (const [k, p] of Object.entries(params)) {
    ret[k] = { name: 'on', params: PD.getDefaultValues((p.map('on') as PD.Group<any>).params) };
  }
  return ret;
}

const registries = (() => {
  const colorThemeRegistry = ColorTheme.createRegistry();
  for (const p of Object.values(BuiltInColorThemes)) colorThemeRegistry.add(p as ColorTheme.Provider<any, any, any>);
  colorThemeRegistry.add(InteractionTypeColorThemeProvider as any);
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

const ctx = { runtime: SyncRuntimeContext, assetManager: new AssetManager() };

async function createRepr(structure: Structure, graphics: any, props: Record<string, any> = {}, sizeTheme = 'uniform') {
  const provider = InteractionsRepresentationProvider;
  const themeProps = {
    colorTheme: { name: 'interaction-type', params: {} },
    sizeTheme: { name: sizeTheme, params: {} },
  };
  await Theme.ensureDependencies(ctx, registries, { structure }, themeProps);
  const repr = provider.factory({ graphics, ...registries }, provider.getParams) as Representation.Any;
  repr.setTheme(Theme.create(registries, { structure }, themeProps));
  const params = PD.getDefaultValues(provider.getParams(registries, structure));
  await repr.createOrUpdate({ ...params, quality: 'low', ...props }, structure).run();
  return repr;
}

function hash(xs: ArrayLike<number>) {
  return (hashFnv32a(xs) >>> 0).toString(16);
}

function lociRows(loci: Loci): number[] {
  if (Interactions.isLoci(loci)) {
    return loci.elements.flatMap((e) => [e.unitA.id, e.indexA, e.unitB.id, e.indexB]);
  } else if (Bridges.isLoci(loci)) {
    return loci.elements.map((e) => e.bridgeIndex);
  }
  return [-1];
}

function markedIndices(ro: GraphicsRenderObject) {
  const marker = (ro.values as any).tMarker.ref.value.array as Uint8Array;
  const count = (ro.values as any).uGroupCount.ref.value * (ro.values as any).instanceCount.ref.value;
  const ret: number[] = [];
  for (let i = 0; i < count; ++i) if (marker[i] !== 0) ret.push(i);
  return ret;
}

/** Positions of the group's geometry: the centroid of its vertices (mesh) or of its dash midpoints (impostors) */
function groupCentroids(ro: GraphicsRenderObject) {
  const v = ro.values as any;
  const groupCount: number = v.uGroupCount.ref.value;
  const sum = new Float64Array(groupCount * 3);
  const n = new Uint32Array(groupCount);
  const group = v.aGroup.ref.value as Float32Array;
  if (ro.type === 'mesh') {
    const pos = v.aPosition.ref.value as Float32Array;
    const index = v.elements.ref.value as Uint32Array;
    const seen = new Uint8Array(group.length);
    for (let k = 0, kl = v.drawCount.ref.value; k < kl; ++k) {
      const i = index[k];
      if (seen[i]) continue;
      seen[i] = 1;
      const g = group[i];
      sum[g * 3] += pos[i * 3];
      sum[g * 3 + 1] += pos[i * 3 + 1];
      sum[g * 3 + 2] += pos[i * 3 + 2];
      n[g]++;
    }
  } else if (ro.type === 'cylinders') {
    const start = v.aStart.ref.value as Float32Array;
    const end = v.aEnd.ref.value as Float32Array;
    // six vertices per cylinder, four triangles
    for (let i = 0, il = (v.drawCount.ref.value / 12) * 6; i < il; i += 6) {
      const g = group[i];
      for (let k = 0; k < 3; ++k) sum[g * 3 + k] += (start[i * 3 + k] + end[i * 3 + k]) / 2;
      n[g]++;
    }
  } else {
    throw new Error(`unexpected render object type '${ro.type}'`);
  }
  const ret = new Map<number, Vec3>();
  for (let g = 0; g < groupCount; ++g) {
    if (n[g]) ret.set(g, Vec3.create(sum[g * 3] / n[g], sum[g * 3 + 1] / n[g], sum[g * 3 + 2] / n[g]));
  }
  return ret;
}

/**
 * Per render object: group and instance count, the drawn groups, the loci each drawn group picks, and the groups each
 * picked loci marks; then the groups a residue loci marks across all render objects
 */
function signature(repr: Representation.Any, structure: Structure) {
  const ros = repr.renderObjects;
  const rows: string[] = [];
  for (const ro of ros) {
    const v = ro.values as any;
    const groupCount: number = v.uGroupCount.ref.value;
    const instanceCount: number = v.instanceCount.ref.value;
    const drawn = [...groupCentroids(ro).keys()];
    const picked: number[] = [];
    const marked: number[] = [];
    for (const instanceId of [0, instanceCount - 1]) {
      for (const groupId of drawn) {
        const loci = repr.getLoci({ objectId: ro.id, instanceId, groupId });
        picked.push(...lociRows(loci));
        repr.mark(loci, MarkerAction.Highlight);
        for (const r of ros) for (const i of markedIndices(r)) marked.push(r.id === ro.id ? i : -1 - i);
        repr.mark(loci, MarkerAction.RemoveHighlight);
      }
    }
    rows.push(
      `groups=${groupCount} instances=${instanceCount} drawn=${drawn.length}:${hash(drawn)} picked=${hash(picked)} marked=${hash(marked)}`,
    );
  }
  rows.sort();

  // residue loci: every residue of the first unit, marked with the 'select' semantics (both features inside)
  const unit = structure.units[0];
  const marked: number[] = [];
  const { residueAtomSegments } = unit.model.atomicHierarchy;
  const { elements } = unit;
  let start = 0;
  for (let i = 1; i <= elements.length; ++i) {
    if (i < elements.length && residueAtomSegments.index[elements[i]] === residueAtomSegments.index[elements[start]])
      continue;
    const indices = new Int32Array(i - start).map((_, k) => start + k);
    const loci = StructureElement.Loci(structure, [{ unit, indices: indices as any }]);
    repr.mark(loci, MarkerAction.Highlight);
    for (const ro of ros) marked.push(ro.id - ros[0].id, ...markedIndices(ro));
    repr.mark(loci, MarkerAction.RemoveHighlight);
    start = i;
  }
  rows.push(`residues=${hash(marked)}`);
  return rows;
}

describe('interactions representation', async () => {
  const m1cbs = await loadCif('1cbs_updated.cif');
  const structure = await StructureSymmetry.buildSymmetryRange(
    Structure.ofModel(m1cbs),
    Vec3.create(0, 0, 0),
    Vec3.create(1, 1, 0),
  ).run();
  await InteractionsProvider.attach(
    ctx,
    structure,
    {
      providers: allOn(ContactProviderParams) as any,
      bridges: allOn(BridgeProviderParams) as any,
    },
    true,
  );

  const expected = [
    'groups=5308 instances=1 drawn=3714:b9d4f0a8 picked=d0bf6c41 marked=b45c7705',
    'groups=614 instances=16 drawn=430:99358b04 picked=14680205 marked=be41d1f9',
    'groups=6628 instances=1 drawn=3470:c0e49c76 picked=e2375f35 marked=b81ec6f5',
    'residues=fe1b5bba',
  ];

  it('mesh visuals: drawn groups, picking and marking', async () => {
    const repr = await createRepr(structure, undefined);
    expect(repr.renderObjects.every((ro) => ro.type === 'mesh')).toBe(true);
    expect(signature(repr, structure)).toEqual(expected);
  });

  it('impostor visuals: the same groups, picking and marking', async () => {
    const repr = await createRepr(structure, graphicsStub);
    expect(repr.renderObjects.every((ro) => ro.type === 'cylinders')).toBe(true);
    expect(signature(repr, structure)).toEqual(expected);
  });

  it('impostor dashes sit where the mesh dashes are, with the same radius', async () => {
    const mesh = await createRepr(structure, undefined);
    const impostors = await createRepr(structure, graphicsStub);
    const { sizeFactor, dashScale, dashCount } = PD.getDefaultValues(
      InteractionsRepresentationProvider.getParams(registries, structure),
    );
    expect(mesh.renderObjects.length).toBe(impostors.renderObjects.length);
    for (let i = 0; i < mesh.renderObjects.length; ++i) {
      const m = groupCentroids(mesh.renderObjects[i]);
      const c = groupCentroids(impostors.renderObjects[i]);
      expect([...c.keys()]).toEqual([...m.keys()]);
      for (const [g, p] of m) expect(Vec3.distance(p, c.get(g)!)).toBeLessThan(1e-3);

      const v = impostors.renderObjects[i].values as any;
      const groups = v.aGroup.ref.value as Float32Array;
      const cylinderCount = v.drawCount.ref.value / 12;
      // one impostor per dash, ceil(dashCount / 2) dashes per drawn half
      expect(cylinderCount).toBe(c.size * Math.ceil(dashCount / 2));
      // buffers allocated once at their final size
      expect(groups.length).toBe(cylinderCount * 6);
      // radius: uniform size 1 × sizeFactor × dashScale, as the mesh bakes it
      expect(v.uSize.ref.value * v.uSizeFactor.ref.value).toBeCloseTo(sizeFactor);
      const scales = (v.aScale.ref.value as Float32Array).subarray(0, cylinderCount * 6);
      expect(scales.every((x) => Math.abs(x - dashScale) < 1e-6)).toBe(true);
    }
  });

  it('switches between impostors and mesh with tryUseImpostor', async () => {
    const repr = await createRepr(structure, graphicsStub, { tryUseImpostor: false });
    expect(repr.renderObjects.every((ro) => ro.type === 'mesh')).toBe(true);
    await repr.createOrUpdate({ tryUseImpostor: true }, structure).run();
    expect(repr.renderObjects.every((ro) => ro.type === 'cylinders')).toBe(true);
    expect(signature(repr, structure)).toEqual(expected);
  });

  it("draws colorMode 'interpolate' as 'default' (the iterators give one location per group)", async () => {
    const repr = await createRepr(structure, graphicsStub, { colorMode: 'interpolate' });
    for (const ro of repr.renderObjects) {
      const v = ro.values as any;
      const modes = (v.aColorMode.ref.value as Float32Array).subarray(0, (v.drawCount.ref.value / 12) * 6);
      expect(modes.every((x) => x === 2)).toBe(true);
    }
  });

  it.each(['physical', 'uncertainty'])(
    "%s size theme: an interaction's size is its smaller end's, as the mesh bakes it",
    async (name) => {
      const theme = registries.sizeThemeRegistry.create(name, { structure }, {});
      const interactions = InteractionsProvider.get(structure).value!;
      const { contacts, unitsFeatures } = interactions;
      const vdw = (unitId: number, feature: number) => {
        const unit = structure.unitMap.get(unitId);
        const f = unitsFeatures.get(unitId);
        const l = StructureElement.Location.create(structure, unit, unit.elements[f.members[f.offsets[feature]]]);
        return theme.size(l);
      };
      const location = Interactions.Location(interactions, structure);
      let different = 0;
      for (const e of contacts.edges) {
        location.element.unitA = structure.unitMap.get(e.unitA);
        location.element.indexA = e.indexA;
        location.element.unitB = structure.unitMap.get(e.unitB);
        location.element.indexB = e.indexB;
        const a = vdw(e.unitA, e.indexA);
        const b = vdw(e.unitB, e.indexB);
        if (a !== b) different++;
        expect(theme.size(location)).toBe(Math.min(a, b));
      }
      expect(different).toBeGreaterThan(0);

      const { bridges } = interactions;
      const bridge = Bridges.Location({ structure, bridges, unitsFeatures });
      for (let i = 0; i < bridges.length; ++i) {
        const b = bridges[i];
        bridge.element.bridgeIndex = i;
        const min = Math.min(
          vdw(b.unitA, b.indexA),
          vdw(b.unitM, b.indexMA),
          vdw(b.unitM, b.indexMB),
          vdw(b.unitB, b.indexB),
        );
        expect(theme.size(bridge)).toBe(min);
      }
    },
  );

  it('updates group sizes when the interactions are recomputed on the same structure', async () => {
    const cell = await StructureSymmetry.buildSymmetryRange(
      Structure.ofModel(m1cbs),
      Vec3.create(0, 0, 0),
      Vec3.create(1, 0, 0),
    ).run();
    await InteractionsProvider.attach(ctx, cell, {}, true);
    const repr = await createRepr(cell, graphicsStub, {}, 'physical');
    await InteractionsProvider.attach(ctx, cell, { providers: allOn(ContactProviderParams) as any }, true);
    await repr.createOrUpdate({}, cell).run();
    const fresh = await createRepr(cell, graphicsStub, {}, 'physical');

    expect(repr.renderObjects.length).toBe(fresh.renderObjects.length);
    for (let i = 0; i < repr.renderObjects.length; ++i) {
      const a = repr.renderObjects[i].values as any;
      const b = fresh.renderObjects[i].values as any;
      const groupCount = b.uGroupCount.ref.value;
      expect(a.uGroupCount.ref.value).toBe(groupCount);
      const sizeA = (a.tSize.ref.value.array as Uint8Array).subarray(0, groupCount * 4);
      const sizeB = (b.tSize.ref.value.array as Uint8Array).subarray(0, groupCount * 4);
      expect(sizeA).toEqual(sizeB);
    }
  });

  it('includeParent: impostors draw, pick and mark as the mesh for every parentDisplay', async () => {
    // a child of the first unit's first half and the whole second unit
    const [u0, u1] = structure.units;
    const half = u0.elements.subarray(0, u0.elements.length >> 1) as any;
    const child = Structure.create([u0.getChild(half), u1], { parent: structure });
    for (const parentDisplay of ['stub', 'full', 'between']) {
      const props = { includeParent: true, parentDisplay };
      const mesh = await createRepr(child, undefined, props);
      const impostors = await createRepr(child, graphicsStub, props);
      expect(impostors.renderObjects.every((ro) => ro.type === 'cylinders')).toBe(true);
      const parent = child.asParent();
      expect(signature(impostors, parent)).toEqual(signature(mesh, parent));
    }
  });
});
