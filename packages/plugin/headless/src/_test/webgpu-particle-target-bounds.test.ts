/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { readFileSync } from 'node:fs';
import { describe, expect, it } from 'vitest';
import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import { Asset } from '@molstar/core/util/assets';
import { getMaxSize } from '@molstar/graphics/geo/geometry/size-data';
import type { GraphicsRenderObject } from '@molstar/graphics/gl/render-object';
import { HeadlessPluginContext } from '@molstar/plugin-headless/context';
import { DefaultPluginSpec } from '@molstar/plugin/default-spec';
import { ParticlesRegistry } from '@molstar/plugin/particles-registry';
import { RelionStarParticlesProvider } from '@molstar/plugin/state/formats/particles/star';
import { ParticleListWithTargets } from '@molstar/plugin/state/transforms/particles/ops';
import { ParticlesRepresentation3D } from '@molstar/plugin/state/transforms/particles/representation';
import { WebGPU } from './webgpu-module.js';

const readExample = (name: string) => readFileSync(new URL(`../../../../../data/examples/${name}`, import.meta.url));

// a structure target's spheres are sized by its elements' radii, not by the representation's per-particle size theme:
// the padding of the bounds (and the sphere raster's size test) follow the elements' sizes
describe.skipIf(!WebGPU.module)('particle structure targets on WebGPU', () => {
  it("pad the bounds by the elements' sizes", async () => {
    const spec = DefaultPluginSpec();
    const plugin = await HeadlessPluginContext.create(
      { webgpu: WebGPU.module! },
      { ...spec, registry: [...(spec.registry ?? []), ...ParticlesRegistry] },
      { width: 32, height: 32 },
      { backend: 'webgpu' },
    );
    try {
      await plugin.init();
      const star = await plugin.builders.data.rawData({
        data: readExample('particles/relion-subtomo.star').toString(),
        label: 'star',
      });
      const { list } = await RelionStarParticlesProvider.parse(plugin, star);
      const targetId = list.cell!.obj!.data.targetInfo.keys().next().value!;
      const file = Asset.File(new File([readExample('1crn.cif')], '1crn.cif'));
      const repr = await plugin
        .build()
        .to(list)
        .apply(ParticleListWithTargets, {
          ...PD.getDefaultValues(ParticleListWithTargets.definition.params!(undefined, plugin)),
          targets: [{ targetId, source: { name: 'file', params: { file, format: 'auto' } } }],
        } as any)
        .apply(ParticlesRepresentation3D, {
          type: { name: 'target', params: {} },
          // the particle's radius, as the Viewer's default visuals size the targets
          sizeTheme: { name: 'particle-size', params: { scale: 1, defaultSize: 50 } },
        } as any)
        .commit();

      const objects = (repr.cell!.obj!.data.repr.renderObjects as GraphicsRenderObject[]).flatMap((o) =>
        'members' in o ? o.members : [o],
      );
      const spheres = objects.filter((o) => o.type === 'spheres');
      expect(spheres.length).toBeGreaterThan(0);
      for (const o of spheres) {
        const values = o.values as any;
        const size = getMaxSize(values);
        // the elements' radii
        expect(size).toBeLessThan(3);
        expect(values.padding.ref.value).toBeCloseTo(size * values.uSizeFactor.ref.value, 5);
      }
    } finally {
      plugin.dispose();
    }
  });
});
