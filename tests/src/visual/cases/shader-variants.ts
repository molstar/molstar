/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import type { GraphicsContext } from '@molstar/graphics/gl/context';
import { isWebGPUContext } from '@molstar/graphics/gl/webgpu/context';
import { MolViewSpec } from '@molstar/mvs/behavior';
import { MVSRuntimeRegistry } from '@molstar/mvs/registry';
import { PluginContext } from '@molstar/plugin/context';
import { DefaultPluginSpec } from '@molstar/plugin/default-spec';
import { PluginSpec } from '@molstar/plugin/spec';
import {
  buildSources,
  type CompileRecord,
  compileSource,
  prepareAllKernels,
  type ReadFile,
  recordCompiles,
  type SourceName,
  SourceNames,
  walkVariants,
} from '../../shader-variants/matrix.js';
import type { BackendProfile } from '../framework/backend.js';
import { type VisualTestCase } from '../framework/test-case.js';

// The variants scenes reach (`.v6/plans/gl/shader-analysis.md`, step 3) compiled on each profile: GLSL by ANGLE and the
// driver, WGSL by Tint. A case per source (`--filter 'shaders/<source>'` runs one: a problem in isolation), one for the
// kernels, and one for the coverage of every program (all sources, minutes per profile)

/** The repository's files, as the visual suite's server serves them from its root */
const readServedFile: ReadFile = async (path) => {
  const response = await fetch(`/${path}`);
  if (!response.ok) throw new Error(`'${path}': ${response.status}`);
  return response.text();
};

/**
 * Whether every program can run on `graphics` (WebGPU with blending into float32 targets, `floatBlend`; WebGL2 with float render
 * targets, several render targets and fragment depth): the profiles where every program has to be compiled; the
 * others lack some passes (OIT, illumination, Hi-Z, GPU surfaces, impostors)
 */
function hasEveryCapability(graphics: GraphicsContext) {
  const { floatRenderTarget, multipleRenderTargets, fragDepth } = graphics.capabilities;
  if (isWebGPUContext(graphics)) return graphics.capabilities.floatBlend;
  return (
    (graphics as { isWebGL2?: boolean }).isWebGL2 === true && floatRenderTarget && multipleRenderTargets && fragDepth
  );
}

/**
 * Runs `compile` with a plugin on a new context of the profile (`createContext`: nothing compiled with it before, so
 * every program is recorded as it is compiled): the compilers' problems, and what `report` adds; a failure that stops
 * `compile` is reported with them
 */
async function withVariantPlugin(
  profile: BackendProfile,
  compile: (plugin: PluginContext, record: CompileRecord) => Promise<void>,
  report: (graphics: GraphicsContext, record: CompileRecord) => string[] = () => [],
) {
  const container = document.createElement('div');
  container.style.cssText = 'position: absolute; left: -10000px; width: 256px; height: 256px';
  const canvas = document.createElement('canvas');
  canvas.width = canvas.height = 256;
  container.appendChild(canvas);
  document.body.appendChild(container);
  const spec = DefaultPluginSpec();
  const plugin = new PluginContext({
    ...spec,
    registry: [...(spec.registry ?? []), MVSRuntimeRegistry],
    behaviors: [...spec.behaviors, PluginSpec.Behavior(MolViewSpec)],
  });
  let record: CompileRecord | undefined;
  try {
    await plugin.init();
    const context = await profile.createContext(canvas, plugin.managers.asset);
    if (!context) return [`no context of profile '${profile.id}'`];
    const { graphics } = context;
    // before the plugin's canvas makes its passes
    record = recordCompiles(graphics);
    await plugin.initViewerAsync(canvas, container, context);
    let failed: string | undefined;
    try {
      await compile(plugin, record);
    } catch (e) {
      failed = `failed: ${e instanceof Error ? e.message : e}`;
    }
    return [
      ...(failed ? [failed] : report(graphics, record)),
      ...(await record.problems()),
      ...((graphics as { takeErrors?: () => string[] }).takeErrors?.() ?? []),
    ];
  } finally {
    record?.stop();
    plugin.dispose();
    container.remove();
  }
}

const sourceCase = (name: SourceName): VisualTestCase => ({
  id: `shaders/${name.replace(' ', '-')}`,
  description: `The variants of the ${name} compile on the profile.`,
  isolate: true,
  scene: async () => [],
  check: (_, profile) =>
    withVariantPlugin(profile, async (plugin, record) => {
      const sources = await buildSources(plugin, readServedFile);
      await compileSource(plugin, sources.find((s) => s.name === name)!, record.count);
    }),
});

/** Cases run only when a filter or tags select them: minutes per profile */
export const ShaderVariantCases: VisualTestCase[] = [
  ...SourceNames.map(sourceCase),
  {
    id: 'shaders/kernels',
    description: 'Every compute kernel compiles on the profile (WebGPU).',
    scene: async () => [],
    check: (_, profile) => withVariantPlugin(profile, (plugin) => prepareAllKernels(plugin.canvas3d!.graphics)),
  },
  {
    id: 'shaders/coverage',
    description: 'Every program is compiled by the variants of all sources on profiles with every capability.',
    scene: async () => [],
    check: (_, profile) =>
      withVariantPlugin(
        profile,
        async (plugin, record) => {
          await walkVariants(plugin, readServedFile, record.count);
          await prepareAllKernels(plugin.canvas3d!.graphics);
        },
        (g, record) =>
          hasEveryCapability(g)
            ? record
                .uncovered()
                // the visual suite's own kernels (`cases/compute.ts`)
                .filter((name) => !name.startsWith('test-'))
                .map((name) => `not compiled: ${name}`)
            : [],
      ),
  },
];
