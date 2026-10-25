/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { readFileSync } from 'node:fs';
import { describe, expect, it } from 'vitest';
import { getSignatureKey, type ObjectSignature } from '@molstar/graphics/gl/precompile';
import { HeadlessPluginContext } from '@molstar/plugin-headless/context';
import { DefaultPluginSpec } from '@molstar/plugin/default-spec';
import { getTreePrecompileList } from '@molstar/plugin/state/precompile';
import { OverpaintStructureRepresentation3DFromScript } from '@molstar/plugin/state/transforms/structure/effects/overpaint';
import { TransparencyStructureRepresentation3DFromScript } from '@molstar/plugin/state/transforms/structure/effects/transparency';
import { WebGPU } from './webgpu-module.js';

const readExample = (name: string) =>
  readFileSync(new URL(`../../../../../data/examples/${name}`, import.meta.url), 'utf8');

/** The ways `r` drew (opaque, transparent, emissive) that its derived signature `d` doesn't allow */
function getMissingFlags(r: ObjectSignature, d: ObjectSignature) {
  return (['opaque', 'transparent', 'emissive'] as const).filter((f) => r[f] && !d[f]);
}

/** The defines (and other fields) in which `r` differs from the closest of `derived`, to name what a derivation misses */
function describeMiss(r: ObjectSignature, derived: readonly ObjectSignature[]) {
  let best: string[] | undefined;
  for (const d of derived) {
    if (d.type !== r.type) continue;
    const diff = Object.keys({ ...r.defines, ...d.defines })
      .filter((k) => r.defines[k] !== d.defines[k])
      .map((k) => `${k}: ${r.defines[k]} (derived ${d.defines[k]})`);
    if (r.doubleSided !== d.doubleSided) diff.push(`doubleSided: ${r.doubleSided}`);
    if (r.merged !== d.merged) diff.push(`merged: ${r.merged}`);
    if (r.linearTextures.join() !== d.linearTextures.join()) diff.push(`linear: ${r.linearTextures}`);
    if (!best || diff.length < best.length) best = diff;
  }
  return `${r.type}${r.merged ? ' merged' : ''}: ${best ? best.join('; ') : 'no derived signature of the kind'}`;
}

const createPlugin = async () => {
  const plugin = await HeadlessPluginContext.create(
    { webgpu: WebGPU.module! },
    DefaultPluginSpec(),
    { width: 32, height: 32 },
    { backend: 'webgpu' },
  );
  await plugin.init();
  return plugin;
};

/** Every structure representation type, with props, themes and layers, on a protein with a ligand and on DNA */
async function addScene(plugin: HeadlessPluginContext) {
  const { representation } = plugin.builders.structure;
  const loadStructure = async (name: string) => {
    const data = await plugin.builders.data.rawData({ data: readExample(name), label: name });
    const trajectory = await plugin.builders.structure.parseTrajectory(data, 'mmcif');
    const model = await plugin.builders.structure.createModel(trajectory);
    return plugin.builders.structure.createStructure(model);
  };

  // every type with its defaults (text needs a canvas module, not available here)
  const protein = await loadStructure('1cbs_updated.cif');
  for (const [type] of plugin.representation.structure.registry.types) {
    if (type === 'label') continue;
    await representation.addRepresentation(protein, { type: type as any });
  }
  // props and layers
  await representation.addRepresentation(protein, { type: 'ball-and-stick', typeParams: { alpha: 0.5 } });
  await representation.addRepresentation(protein, {
    type: 'molecular-surface',
    typeParams: { doubleSided: true, transparentBackfaces: 'on', quality: 'low' },
    color: 'secondary-structure',
  });
  await representation.addRepresentation(protein, {
    type: 'gaussian-surface',
    typeParams: { xrayShaded: true },
    color: 'uniform',
  });
  await representation.addRepresentation(protein, { type: 'cartoon', color: 'sequence-id' });
  const spacefill = await representation.addRepresentation(protein, { type: 'spacefill', color: 'chain-id' });
  await plugin
    .build()
    .to(spacefill)
    .apply(OverpaintStructureRepresentation3DFromScript)
    .to(spacefill)
    .apply(TransparencyStructureRepresentation3DFromScript)
    .commit();

  // nucleic acids
  const dna = await loadStructure('1bna_confal_pyramids.cif');
  for (const type of ['cartoon', 'ball-and-stick', 'line'] as const) {
    await representation.addRepresentation(dna, { type });
  }

  // merged (`renderMerged` on: 'auto' is derived as unmerged); one chain in 40 copies, merged alone to be culled
  // on the GPU
  const loadPdb = async (pdb: string, assembly: boolean) => {
    const data = await plugin.builders.data.rawData({ data: pdb, label: 'pdb' });
    const trajectory = await plugin.builders.structure.parseTrajectory(data, 'pdb');
    const model = await plugin.builders.structure.createModel(trajectory);
    return plugin.builders.structure.createStructure(
      model,
      assembly ? { name: 'assembly', params: { id: '1' } } : undefined,
    );
  };
  const chains = await loadPdb(chainsPdb(32, 0), false);
  await representation.addRepresentation(chains, { type: 'spacefill', typeParams: { renderMerged: true } });
  await representation.addRepresentation(chains, { type: 'ball-and-stick', typeParams: { renderMerged: true } });
  const copies = await loadPdb(copiesPdb(40), true);
  await representation.addRepresentation(copies, { type: 'spacefill', typeParams: { gpuCulling: true } });
  plugin.canvas3d!.commit(true);
}

const f = (v: number) => v.toFixed(3).padStart(8);
/** `count` chains of four alanines' CA atoms, from `copies` assembly operators (none: 0) */
function chainsPdb(count: number, copies: number) {
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
  const biomt: string[] = [];
  if (copies) {
    biomt.push('REMARK 350 BIOMOLECULE: 1', 'REMARK 350 APPLY THE FOLLOWING TO CHAINS: A');
    for (let k = 0; k < copies; ++k) {
      const n = String(k + 1).padStart(4);
      biomt.push(`REMARK 350   BIOMT1 ${n}  1.000000  0.000000  0.000000 ${f((k % 8) * 20)}`);
      biomt.push(`REMARK 350   BIOMT2 ${n}  0.000000  1.000000  0.000000 ${f(Math.floor(k / 8) * 20)}`);
      biomt.push(`REMARK 350   BIOMT3 ${n}  0.000000  0.000000  1.000000 ${f(0)}`);
    }
  }
  return `${[...biomt, ...atoms].join('\n')}\nEND\n`;
}
const copiesPdb = (copies: number) => chainsPdb(1, copies);

describe.skipIf(!WebGPU.module)('precompile lists from params on WebGPU', () => {
  it('cover what the representations draw, and prepare it before the objects exist', async () => {
    const a = await createPlugin();
    const b = await createPlugin();
    try {
      const canvas = a.canvas3d!;
      canvas.resetPrecompileList();
      await addScene(a);

      const recorded = canvas.getPrecompileList().objects;
      const derived = getTreePrecompileList(a, a.state.data.tree).objects;
      const derivedKeys = new Set(derived.map(getSignatureKey));
      const missing = recorded.filter((r) => !derivedKeys.has(getSignatureKey(r))).map((r) => describeMiss(r, derived));
      expect(recorded.length, `recorded ${recorded.length}, derived ${derived.length}`).toBeGreaterThan(10);
      expect(missing).toEqual([]);
      const byKey = new Map(derived.map((d) => [getSignatureKey(d), d]));
      const flags = recorded
        .map((r) => [r.type, getMissingFlags(r, byKey.get(getSignatureKey(r))!)] as const)
        .filter(([, f]) => f.length);
      expect(flags).toEqual([]);
      // merged objects, one of them culled on the GPU: derived as well, the GPU culling's kernels prepared with them
      expect(recorded.some((r) => r.merged)).toBe(true);
      const gpuCulled = recorded.filter((r) => r.gpuCulled);
      expect(gpuCulled.length).toBeGreaterThan(0);
      expect(gpuCulled.filter((r) => !byKey.get(getSignatureKey(r))?.gpuCulled).map((r) => r.type)).toEqual([]);
      // the over-approximation stays small: a few derived per recorded signature
      expect(derived.length, `recorded ${recorded.length}, derived ${derived.length}`).toBeLessThan(
        recorded.length * 4,
      );

      // in a new context, the derived list prepared first: the same scene, its frames and first interaction compile no
      // render pipeline more
      const pending = b.canvas3d!.graphics.resources.pending;
      await b.canvas3d!.prepare({ interaction: true });
      await b.canvas3d!.precompile({ objects: derived });
      const before = pending.trackedCount;
      await addScene(b);
      await b.canvas3d!.prepare({ interaction: true });
      await pending.whenReady();
      // render pipelines only: the GPU surfaces' compute kernels and readbacks are their own (prepared as they compute)
      expect(pending.trackedSince(before).filter((label) => label.startsWith('pipeline '))).toEqual([]);
    } finally {
      a.dispose();
      b.dispose();
    }
  }, 120_000);
});
