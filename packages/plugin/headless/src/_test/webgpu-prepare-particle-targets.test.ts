/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { readFileSync } from 'node:fs';
import { describe, expect, it } from 'vitest';
import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import { Asset } from '@molstar/core/util/assets';
import { now } from '@molstar/core/util/now';
import { getSignatureKey, type ObjectSignature, type PrecompileList } from '@molstar/graphics/gl/precompile';
import { HeadlessPluginContext } from '@molstar/plugin-headless/context';
import { DefaultPluginSpec } from '@molstar/plugin/default-spec';
import { ParticlesRegistry } from '@molstar/plugin/particles-registry';
import { RelionStarParticlesProvider } from '@molstar/plugin/state/formats/particles/star';
import { ParticleListWithTargets } from '@molstar/plugin/state/transforms/particles/ops';
import { ParticlesRepresentation3D } from '@molstar/plugin/state/transforms/particles/representation';
import { WebGPU } from './webgpu-module.js';

const readExample = (name: string) => readFileSync(new URL(`../../../../../data/examples/${name}`, import.meta.url));

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

/** A unit cube as OBJ */
const CubeObj = [
  'v 0 0 0',
  'v 1 0 0',
  'v 1 1 0',
  'v 0 1 0',
  'v 0 0 1',
  'v 1 0 1',
  'v 1 1 1',
  'v 0 1 1',
  'f 1 3 2',
  'f 1 4 3',
  'f 5 6 7',
  'f 5 7 8',
  'f 1 2 6',
  'f 1 6 5',
  'f 2 3 7',
  'f 2 7 6',
  'f 3 4 8',
  'f 3 8 7',
  'f 4 1 5',
  'f 4 5 8',
  '',
].join('\n');

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
    if (r.linearTextures.join() !== d.linearTextures.join()) diff.push(`linear: ${r.linearTextures}`);
    if (!best || diff.length < best.length) best = diff;
  }
  return `${r.type}${r.merged ? ' merged' : ''}: ${best ? best.join('; ') : 'no derived signature of the kind'}`;
}

const Targets = [
  { kind: 'structure', file: () => new File([readExample('1crn.cif')], '1crn.cif') },
  { kind: 'volume', file: () => new File([ballDx(12)], 'ball.dx') },
  { kind: 'shape', file: () => new File([CubeObj], 'cube.obj') },
];

async function createPlugin() {
  const spec = DefaultPluginSpec();
  const plugin = await HeadlessPluginContext.create(
    { webgpu: WebGPU.module! },
    { ...spec, registry: [...(spec.registry ?? []), ...ParticlesRegistry] },
    { width: 32, height: 32 },
    { backend: 'webgpu' },
  );
  await plugin.init();
  return plugin;
}

/** RELION particles whose first target is `file`'s object, drawn as the target representation */
async function addTargets(plugin: HeadlessPluginContext, file: File) {
  const star = await plugin.builders.data.rawData({
    data: readExample('particles/relion-subtomo.star').toString(),
    label: 'star',
  });
  const { list } = await RelionStarParticlesProvider.parse(plugin, star);
  const particles = list.cell!.obj!.data;
  const targetId = particles.targetInfo.keys().next().value!;
  await plugin
    .build()
    .to(list)
    .apply(ParticleListWithTargets, {
      ...PD.getDefaultValues(ParticleListWithTargets.definition.params!(undefined, plugin)),
      targets: [{ targetId, source: { name: 'file', params: { file: Asset.File(file), format: 'auto' } } }],
    } as any)
    .apply(ParticlesRepresentation3D, { type: { name: 'target', params: {} } } as any)
    .commit();
}

// the particle target representation derives its signatures from its particles' targets (a structure, a volume, a
// shape each), so its transform prepares them while it builds (`prepareRepresentation`)
describe.skipIf(!WebGPU.module)('particle targets prepared while they build on WebGPU', () => {
  for (const target of Targets) {
    it(`cover what a ${target.kind} target draws, and leave its first frame nothing to compile`, async () => {
      const plugin = await createPlugin();
      try {
        const canvas = plugin.canvas3d!;
        const prepared: PrecompileList[] = [];
        const prepareObjects = canvas.prepareObjects.bind(canvas);
        canvas.prepareObjects = (owner, list) => {
          prepared.push(list);
          prepareObjects(owner, list);
        };
        const { pending } = canvas.graphics.resources;
        let time = now();
        const frame = () => {
          canvas.requestDraw();
          canvas.tick((time += 1000) as now.Timestamp, { isSynchronous: true });
        };
        frame();
        await pending.whenReady();

        canvas.resetPrecompileList();
        await addTargets(plugin, target.file());
        await pending.whenReady();
        const before = pending.trackedCount;
        frame();
        const compiled = pending.trackedSince(before).filter((label) => label.startsWith('pipeline '));

        const recorded = canvas.getPrecompileList().objects;
        const derived = prepared.flatMap((l) => l.objects);
        const keys = new Set(derived.map(getSignatureKey));
        expect(recorded.length).toBeGreaterThan(0);
        expect(recorded.filter((r) => !keys.has(getSignatureKey(r))).map((r) => describeMiss(r, derived))).toEqual([]);
        expect(compiled).toEqual([]);
        await pending.whenReady();
      } finally {
        plugin.dispose();
      }
    }, 120_000);
  }
});
