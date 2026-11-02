/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { readFileSync } from 'node:fs';
import { describe, expect, it } from 'vitest';
import { now } from '@molstar/core/util/now';
import { getSignatureKey, type ObjectSignature, type PrecompileList } from '@molstar/graphics/gl/precompile';
import { HeadlessPluginContext } from '@molstar/plugin-headless/context';
import { SpacefillParticlesRepresentationProvider } from '@molstar/graphics/repr/particles/representation/spacefill';
import { DefaultPluginSpec } from '@molstar/plugin/default-spec';
import { ParticlesRegistry } from '@molstar/plugin/particles-registry';
import { RelionStarParticlesProvider } from '@molstar/plugin/state/formats/particles/star';
import { DxProvider } from '@molstar/plugin/state/formats/volume/dx';
import { getTreePrecompileList } from '@molstar/plugin/state/precompile';
import { ParticlesRepresentation3D } from '@molstar/plugin/state/transforms/particles/representation';
import { WebGPU } from './webgpu-module.js';

const readExample = (name: string) =>
  readFileSync(new URL(`../../../../../data/examples/${name}`, import.meta.url), 'utf8');

const createPlugin = async () => {
  const spec = DefaultPluginSpec();
  const plugin = await HeadlessPluginContext.create(
    { webgpu: WebGPU.module! },
    { ...spec, registry: [...(spec.registry ?? []), ...ParticlesRegistry] },
    { width: 32, height: 32 },
    { backend: 'webgpu' },
  );
  await plugin.init();
  return plugin;
};

/** What the representation transforms prepare while they build (`Canvas3D.prepareObjects`), per call */
function capturePrepared(plugin: HeadlessPluginContext, prepare = true) {
  const canvas = plugin.canvas3d!;
  const lists: PrecompileList[] = [];
  const prepareObjects = canvas.prepareObjects.bind(canvas);
  canvas.prepareObjects = (owner, list) => {
    lists.push(list);
    if (prepare) prepareObjects(owner, list);
  };
  return lists;
}

async function loadStructure(plugin: HeadlessPluginContext, name: string) {
  const data = await plugin.builders.data.rawData({ data: readExample(name), label: name });
  const trajectory = await plugin.builders.structure.parseTrajectory(data, 'mmcif');
  const model = await plugin.builders.structure.createModel(trajectory);
  return plugin.builders.structure.createStructure(model);
}

/** An OpenDX grid of `n`³ points: a ball's density */
function ballDx(n: number) {
  const values: string[] = [];
  const c = (n - 1) / 2;
  for (let i = 0; i < n; ++i) {
    for (let j = 0; j < n; ++j) {
      for (let k = 0; k < n; ++k) values.push((c - Math.hypot(i - c, j - c, k - c)).toFixed(3));
    }
  }
  return [
    `object 1 class gridpositions counts ${n} ${n} ${n}`,
    'origin 0 0 0',
    'delta 1 0 0',
    'delta 0 1 0',
    'delta 0 0 1',
    `object 2 class gridconnections counts ${n} ${n} ${n}`,
    `object 3 class array type double rank 0 items ${n * n * n} data follows`,
    ...values,
    '',
  ].join('\n');
}

/** `count` chains of four alanines' CA atoms */
function chainsPdb(count: number) {
  const f = (v: number) => v.toFixed(3).padStart(8);
  const ids = 'ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz';
  const atoms: string[] = [];
  for (let c = 0; c < count; ++c) {
    for (let i = 0; i < 4; ++i) {
      const serial = String(c * 4 + i + 1).padStart(5);
      atoms.push(
        `ATOM  ${serial}  CA  ALA ${ids[c]}${String(i + 1).padStart(4)}    ${f(i * 3.8)}${f(c * 8)}${f(0)}  1.00  0.00           C`,
      );
    }
  }
  return `${atoms.join('\n')}\nEND\n`;
}

/** The defines in which `r` differs from the closest of `derived`, to name what a derivation misses */
function describeMiss(r: ObjectSignature, derived: readonly ObjectSignature[]) {
  let best: string[] | undefined;
  for (const d of derived) {
    if (d.type !== r.type) continue;
    const diff = Object.keys({ ...r.defines, ...d.defines })
      .filter((k) => r.defines[k] !== d.defines[k])
      .map((k) => `${k}: ${r.defines[k]} (derived ${d.defines[k]})`);
    if (r.doubleSided !== d.doubleSided) diff.push(`doubleSided: ${r.doubleSided}`);
    if (r.merged !== d.merged) diff.push(`merged: ${r.merged}`);
    if (!best || diff.length < best.length) best = diff;
  }
  return `${r.type}${r.merged ? ' merged' : ''}: ${best ? best.join('; ') : 'no derived signature of the kind'}`;
}

// a representation transform derives its signatures with its data and themes before it builds, and the canvas compiles
// them meanwhile: the data resolves what params alone take every way (themes' granularity, 'auto' quality, visuals,
// merging), and the first frame after the commit finds its pipelines
describe.skipIf(!WebGPU.module)('representations prepared while they build on WebGPU', () => {
  it('cover what they draw, also what params alone take otherwise (auto merging)', async () => {
    const plugin = await createPlugin();
    try {
      const prepared = capturePrepared(plugin, false);
      const canvas = plugin.canvas3d!;
      canvas.resetPrecompileList();
      const { representation } = plugin.builders.structure;
      const protein = await loadStructure(plugin, '1cbs_updated.cif');
      // every type with its defaults (text needs a canvas module, not available here), and some props and themes
      for (const [type] of plugin.representation.structure.registry.types) {
        if (type === 'label') continue;
        await representation.addRepresentation(protein, { type: type as any });
      }
      await representation.addRepresentation(protein, { type: 'ball-and-stick', typeParams: { alpha: 0.5 } });
      await representation.addRepresentation(protein, {
        type: 'molecular-surface',
        typeParams: { doubleSided: true, quality: 'low' },
        color: 'secondary-structure',
      });
      await representation.addRepresentation(protein, { type: 'cartoon', color: 'sequence-id' });
      const dna = await loadStructure(plugin, '1bna_confal_pyramids.cif');
      for (const type of ['cartoon', 'ball-and-stick'] as const) await representation.addRepresentation(dna, { type });
      // 32 chains: 'auto' merging merges them (above `AutoRenderMergedThreshold`), which params alone take as unmerged
      const chains = await plugin.builders.data.rawData({ data: chainsPdb(32), label: 'chains' });
      const trajectory = await plugin.builders.structure.parseTrajectory(chains, 'pdb');
      const model = await plugin.builders.structure.createModel(trajectory);
      const structure = await plugin.builders.structure.createStructure(model);
      await representation.addRepresentation(structure, { type: 'spacefill' });
      canvas.commit(true);

      const recorded = canvas.getPrecompileList().objects;
      const withData = prepared.flatMap((l) => l.objects);
      const keys = new Set(withData.map(getSignatureKey));
      const missing = recorded.filter((r) => !keys.has(getSignatureKey(r))).map((r) => describeMiss(r, withData));
      expect(recorded.length).toBeGreaterThan(10);
      expect(missing).toEqual([]);
      expect(recorded.some((r) => r.merged)).toBe(true);
      const fromParams = getTreePrecompileList(plugin, plugin.state.data.tree).objects;
      const paramsKeys = new Set(fromParams.map(getSignatureKey));
      const missingFromParams = recorded.filter((r) => !paramsKeys.has(getSignatureKey(r)));
      expect(missingFromParams.map((r) => `${r.type}${r.merged ? ' merged' : ''}`)).toEqual(['spheres merged']);
    } finally {
      plugin.dispose();
    }
  }, 120_000);

  it('cover what volume and particle representations draw', async () => {
    const plugin = await createPlugin();
    try {
      const prepared = capturePrepared(plugin, false);
      const canvas = plugin.canvas3d!;
      canvas.resetPrecompileList();
      // a volume with the format's visuals (an isosurface)
      const dx = await plugin.builders.data.rawData({ data: ballDx(16), label: 'ball.dx' });
      const volume = await DxProvider.parse(plugin, dx);
      await DxProvider.visuals?.(plugin, volume);
      // particles as spacefill
      const star = await plugin.builders.data.rawData({
        data: readExample('particles/relion-subtomo.star'),
        label: 'star',
      });
      const { list } = await RelionStarParticlesProvider.parse(plugin, star);
      await plugin
        .build()
        .to(list)
        .apply(ParticlesRepresentation3D, {
          type: { name: SpacefillParticlesRepresentationProvider.name, params: {} },
        } as any)
        .commit();
      canvas.commit(true);

      const recorded = canvas.getPrecompileList().objects;
      const withData = prepared.flatMap((l) => l.objects);
      const keys = new Set(withData.map(getSignatureKey));
      const missing = recorded.filter((r) => !keys.has(getSignatureKey(r))).map((r) => describeMiss(r, withData));
      expect([...new Set(recorded.map((r) => r.type))].sort()).toEqual(['spheres', 'texture-mesh']);
      expect(missing).toEqual([]);
    } finally {
      plugin.dispose();
    }
  }, 120_000);

  for (const prepare of [true, false]) {
    it(`${prepare ? 'leave' : 'without it, leave'} the first frame after the commit ${prepare ? 'nothing' : 'pipelines'} to compile`, async () => {
      const plugin = await createPlugin();
      try {
        capturePrepared(plugin, prepare);
        const canvas = plugin.canvas3d!;
        const { pending } = canvas.graphics.resources;
        let time = now();
        const frame = () => {
          canvas.requestDraw();
          canvas.tick((time += 1000) as now.Timestamp, { isSynchronous: true });
        };
        frame();
        await pending.whenReady();
        const protein = await loadStructure(plugin, '1cbs_updated.cif');
        await plugin.builders.structure.representation.addRepresentation(protein, {
          type: 'cartoon',
          color: 'secondary-structure',
        });
        await plugin.builders.structure.representation.addRepresentation(protein, { type: 'ball-and-stick' });
        // what was prepared while they built is ready by now (Dawn compiles one after the other)
        await pending.whenReady();
        const before = pending.trackedCount;
        frame();
        const compiled = pending.trackedSince(before).filter((label) => label.startsWith('pipeline '));
        if (prepare) expect(compiled).toEqual([]);
        else expect(compiled.length).toBeGreaterThan(0);
        await pending.whenReady();
      } finally {
        plugin.dispose();
      }
    }, 120_000);
  }
});
