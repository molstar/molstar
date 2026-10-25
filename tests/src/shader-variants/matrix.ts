/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Color } from '@molstar/core/util/color';
import { ColorNames } from '@molstar/core/util/color/names';
import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import type { Canvas3DPrecompileProps } from '@molstar/graphics/canvas3d/canvas3d';
import { Canvas3DParams } from '@molstar/graphics/canvas3d/canvas3d';
import { Spheres } from '@molstar/graphics/geo/geometry/spheres/spheres';
import type { ObjectSignature, PrecompileList } from '@molstar/graphics/gl/precompile';
import { getRenderableProgramSource } from '@molstar/graphics/gl/render-object';
import type { RenderableSchema } from '@molstar/graphics/gl/renderable/schema';
import { SpacefillRepresentationProvider } from '@molstar/graphics/repr/structure/representation/spacefill';
import { MVSData } from '@molstar/mvs-builder/mvs-data';
import { loadMVS } from '@molstar/mvs/load';
import { now } from '@molstar/core/util/now';
import type { TransparencyMode } from '@molstar/graphics/canvas3d/passes/draw';
import type { GraphicsContext } from '@molstar/graphics/gl/context';
import { getColorSmoothingKernels } from '@molstar/graphics/gl/compute/color-smoothing';
import { GaussianDensityKernels } from '@molstar/graphics/gl/compute/gaussian-density';
import { getIsosurfaceKernels } from '@molstar/graphics/gl/compute/marching-cubes/compute';
import { PrefixSumKernels } from '@molstar/graphics/gl/compute/prefix-sum';
import { MergedCullKernels } from '@molstar/graphics/gl/renderable/merged-culler';
import type { WebGLContext } from '@molstar/graphics/gl/webgl/context';
import { getGlslProgramNames } from '@molstar/graphics/gl/webgl/programs';
import { isWebGPUContext, type WebGPUContext } from '@molstar/graphics/gl/webgpu/context';
import { getWgslKernelNames } from '@molstar/graphics/gl/webgpu/shader/kernels';
import { getWgslProgramNames } from '@molstar/graphics/gl/webgpu/shader/programs';
import type { PluginContext } from '@molstar/plugin/context';
import { createStructureRepresentationParams } from '@molstar/plugin/state/helpers/structure-representation-params';
import { createVolumeRepresentationParams } from '@molstar/plugin/state/helpers/volume-representation-params';
import { getSnapshotsPrecompileList, getTreePrecompileList } from '@molstar/plugin/state/precompile';
import { uniqueSignatures as unique } from '@molstar/graphics/repr/signatures';
import { StructureRepresentation3D } from '@molstar/plugin/state/transforms/structure/representation';
import { VolumeRepresentation3D } from '@molstar/plugin/state/transforms/volume/representation';

// The variants of the shaders that scenes reach (`.v6/plans/gl/shader-analysis.md`, phase 1): the precompile lists
// (object signatures, `getTreePrecompileList`) and the canvas props of three sources, the default presets, the
// MolViewSpec examples and mesoscale scenes, then each define and each canvas option varied alone; and a walk that
// compiles them on a plugin's canvas. For Node (the Dawn and native `gl` tests) and the browser (the visual suite):
// files are read through `ReadFile`.

/** Reads a file of the repository by its path from the root (Node: the file system; the browser: the served files) */
export type ReadFile = (path: string) => Promise<string>;

/** A source of variants: what its scenes draw (object signatures) and the canvas props they draw with */
export interface VariantSource {
  readonly name: string;
  readonly list: PrecompileList;
  readonly props: readonly Canvas3DPrecompileProps[];
}

/** The transparency modes a canvas draws in, each a context of its own for the programs */
export const TransparencyModes = ['blended', 'wboit', 'dpoit'] as const;
const withTransparencyModes = (props: Canvas3DPrecompileProps) =>
  TransparencyModes.map((transparency) => ({ ...props, transparency }));

async function loadStructure(plugin: PluginContext, read: ReadFile, name: string, format: 'mmcif' = 'mmcif') {
  const data = await plugin.builders.data.rawData({ data: await read(`data/examples/${name}`), label: name });
  const trajectory = await plugin.builders.structure.parseTrajectory(data, format);
  const model = await plugin.builders.structure.createModel(trajectory);
  return plugin.builders.structure.createStructure(model);
}

/**
 * The default presets: every structure representation preset that applies, on a protein, a nucleic acid and a coarse
 * model; every structure and volume representation type with its defaults (an uncommitted tree: derived from params,
 * nothing built)
 */
export async function getPresetSource(plugin: PluginContext, read: ReadFile): Promise<VariantSource> {
  const structures = [
    await loadStructure(plugin, read, '1cbs_updated.cif'),
    await loadStructure(plugin, read, '1bna_confal_pyramids.cif'),
    await loadStructure(plugin, read, 'ihm-coarse-8zzh.cif'),
  ];
  const { representation } = plugin.builders.structure;
  for (const structure of structures) {
    for (const preset of representation.getPresets(structure.obj)) {
      await representation.applyPreset(structure.ref, preset);
    }
  }
  const signatures: ObjectSignature[] = [...getTreePrecompileList(plugin, plugin.state.data.tree).objects];

  const builder = plugin.build();
  const [protein] = structures;
  Structures.set(plugin, protein.ref);
  for (const [type] of plugin.representation.structure.registry.types) {
    const params = createStructureRepresentationParams(plugin, protein.data, { type: type as any });
    builder.to(protein.ref).apply(StructureRepresentation3D, params);
  }
  for (const [type] of plugin.representation.volume.registry.types) {
    // under the root: the derivation reads the params (and the parent's), not the data
    builder
      .toRoot()
      .apply(VolumeRepresentation3D, createVolumeRepresentationParams(plugin, undefined, { type: type as any }));
  }
  signatures.push(...getTreePrecompileList(plugin, builder.getTree()).objects);
  return { name: 'presets', list: { objects: unique(signatures) }, props: withTransparencyModes({}) };
}

/**
 * The MolViewSpec documents of the repository (`data/examples/mvs`, the runtime's test fixtures), by path: a list, as
 * the browser can't list directories (the Node test checks it against them)
 */
export const MVSDocuments = [
  'data/examples/mvs/1cbs-focus.mvsj',
  'data/examples/mvs/1cbs.mvsj',
  'data/examples/mvs/1h9t_domain_colors.mvsj',
  'data/examples/mvs/1h9t_domain_labels.mvsj',
  'data/examples/mvs/kinase-story.mvsj',
  'packages/mvs/runtime/src/_test/fixtures/exosome.mvsj',
  'packages/mvs/runtime/src/_test/fixtures/features.mvsj',
  'packages/mvs/runtime/src/_test/fixtures/motm1.mvsj',
] as const;

/**
 * The MolViewSpec examples: each document's snapshots (and the frames of their transitions) and the canvas props they
 * set (`getSnapshotsPrecompileList`). Loaded with every download pointed at a missing local file: the snapshots are
 * made before anything is applied, and nothing reaches the network. Needs the MolViewSpec runtime in the plugin
 */
export async function getMVSSource(plugin: PluginContext, read: ReadFile): Promise<VariantSource> {
  const signatures: ObjectSignature[] = [];
  const props: Canvas3DPrecompileProps[] = [];
  for (const name of MVSDocuments) {
    const document = JSON.parse(await read(name));
    const offline = (node: any) => {
      if (node.kind === 'download') node.params = { ...node.params, url: 'file:///molstar-shader-variants/missing' };
      node.children?.forEach(offline);
    };
    if (document.root) offline(document.root);
    for (const snapshot of document.snapshots ?? []) offline(snapshot.root);
    plugin.managers.snapshot.clear();
    try {
      await loadMVS(plugin, MVSData.fromMVSJ(JSON.stringify(document)), { doNotReportErrors: true, keepCamera: true });
    } catch {
      // applying the first snapshot fails without its data; the snapshots are made by then
    }
    const snapshots = plugin.managers.snapshot.state.entries.map((e) => e.snapshot).toArray();
    if (snapshots.length === 0) throw new Error(`no snapshots from '${name}'`);
    const derived = getSnapshotsPrecompileList(plugin, snapshots);
    signatures.push(...derived.list.objects);
    for (const p of derived.props) props.push(...withTransparencyModes(p));
  }
  return { name: 'mvs', list: { objects: unique(signatures) }, props };
}

/**
 * Mesoscale scenes as the Mesoscale Explorer builds them: its spacefill per graphics mode (`getGraphicsModeProps` and
 * `getSpacefillParams` of `apps/mesoscale-explorer`, from the public presets they use, `Spheres.LodPixelPresets`), with
 * either clip variant, and its canvas props (Hi-Z, no fog, sharpening below full resolution, illumination as an option)
 */
export function getMesoscaleSource(plugin: PluginContext, structureRef: string): VariantSource {
  const builder = plugin.build();
  const modes = ['ultra', 'quality', 'balanced', 'performance'] as const;
  for (const mode of modes) {
    const { lodLevels, lodPixelSize } = Spheres.LodPixelPresets[mode];
    for (const variant of ['instance', 'pixel'] as const) {
      builder.to(structureRef).apply(StructureRepresentation3D, {
        type: {
          name: 'spacefill',
          params: {
            ...SpacefillRepresentationProvider.defaultValues,
            ignoreHydrogens: true,
            instanceGranularity: true,
            ignoreLight: true,
            lodLevels,
            lodPixelSize,
            quality: 'lowest',
            clip: { variant, objects: [] },
            clipPrimitive: true,
            approximate: mode !== 'quality' && mode !== 'ultra',
            alphaThickness: mode === 'performance' ? 15 : 12,
            interior: {
              color: Color.fromNormalizedRgb(0, 0, 0),
              colorStrength: 0.15,
              substance: { metalness: 0, roughness: 1, bumpiness: 0 },
              substanceStrength: 1,
            },
          },
        },
        colorTheme: { name: 'uniform', params: { value: ColorNames.lightgray, lightness: 0, saturation: 0 } },
        sizeTheme: { name: 'physical', params: { scale: 1 } },
      });
    }
  }
  const list = getTreePrecompileList(plugin, builder.getTree());
  const base: Canvas3DPrecompileProps = {
    renderer: { ...PD.getDefaultValues(Canvas3DParams).renderer, backgroundColor: Color(0x101010) },
    cameraFog: { name: 'off', params: {} },
    hiZ: { ...PD.getDefaultValues(Canvas3DParams).hiZ, enabled: true },
  };
  const defaults = PD.getDefaultValues(Canvas3DParams);
  const sharpening = {
    postprocessing: {
      ...defaults.postprocessing,
      sharpening: { name: 'on' as const, params: { sharpness: 0.5, denoise: true } },
    },
  };
  const illumination = { illumination: { ...defaults.illumination, enabled: true } };
  const props = [base, { ...base, ...sharpening }, { ...base, ...illumination }].flatMap(withTransparencyModes);
  return { name: 'mesoscale', list: { objects: unique(list.objects) }, props };
}

/**
 * Each define of each object kind of `sources` varied alone, from the first signature of the kind: booleans flipped,
 * strings through the values the sources reach for the kind (the derivation enumerates every value data decides, e.g.
 * color and layer types; another option is a variant no object of the kind draws, as a direct-volume color type on a
 * mesh), number defines (counts, sizes) as they are; a palette with a color type per group, clip objects with a count
 */
export function getDefineAxisSource(sources: readonly VariantSource[]): VariantSource {
  const baselines = new Map<string, ObjectSignature>();
  // the values each string define takes in the sources, by kind
  const reached = new Map<string, Map<string, Set<string>>>();
  for (const s of sources) {
    for (const o of s.list.objects) {
      const kind = `${o.type}|${o.merged}`;
      if (!baselines.has(kind)) baselines.set(kind, o);
      let values = reached.get(kind);
      if (!values) reached.set(kind, (values = new Map()));
      for (const [name, value] of Object.entries(o.defines)) {
        if (typeof value !== 'string') continue;
        let set = values.get(name);
        if (!set) values.set(name, (set = new Set()));
        set.add(value);
      }
    }
  }
  const signatures: ObjectSignature[] = [];
  for (const [kind, base] of baselines) {
    const schema = getRenderableProgramSource(base)?.schema as RenderableSchema | undefined;
    if (!schema) continue;
    for (const [name, value] of Object.entries(base.defines)) {
      const spec = schema[name];
      if (!spec || spec.type !== 'define') continue;
      const values: (string | boolean)[] =
        spec.kind === 'boolean' ? [!value] : spec.kind === 'string' ? [...(reached.get(kind)?.get(name) ?? [])] : [];
      for (const v of values) {
        if (v === value) continue;
        const defines = { ...base.defines, [name]: v };
        // palette indices are colors per group or instance: no theme gives them as a uniform color
        if (name === 'dUsePalette' && v === true && defines.dColorType === 'uniform') defines.dColorType = 'texture';
        // clip objects and their count go together (`dClipObjects` is a count above zero)
        if (name === 'dClipObjects') defines.dClipObjectCount = v ? Math.max(1, Number(defines.dClipObjectCount)) : 0;
        signatures.push({ ...base, defines });
      }
    }
  }
  // in one transparency mode: the modes select the color variants only, which the sources compile in each mode
  return { name: 'define axes', list: { objects: unique(signatures) }, props: [{ transparency: 'blended' }] };
}

/** A mapped postprocessing param */
const mapped = (name: string) => (Canvas3DParams.postprocessing.params as PD.Params)[name] as PD.Mapped<any>;

/** The canvas groups whose options select programs or their defines */
const CanvasAxisGroups = ['postprocessing', 'multiSample', 'marking', 'illumination', 'hiZ', 'renderer'] as const;

/**
 * Each canvas option varied alone from the defaults: every option of a mapped or select param (with the option's
 * defaults), booleans flipped, in the groups of `CanvasAxisGroups` (nested groups too) and the top-level fog and
 * transparent background; the dimming of what is not marked (a number), illumination with multi-sampling and sharpening
 * with each antialiasing (two options) besides; drawn with the objects of `sources`
 */
export function getCanvasAxisSource(sources: readonly VariantSource[]): VariantSource {
  const defaults = PD.getDefaultValues(Canvas3DParams) as any;
  const props: Canvas3DPrecompileProps[] = [];
  const vary = (params: PD.Params, values: any, wrap: (v: any) => Canvas3DPrecompileProps) => {
    for (const [key, p] of Object.entries(params)) {
      if (p.type === 'mapped') {
        for (const [name] of p.select.options) {
          if (name === values[key].name) continue;
          props.push(wrap({ ...values, [key]: { name, params: p.map(name).defaultValue } }));
        }
      } else if (p.type === 'select') {
        for (const [option] of p.options) {
          if (option !== values[key]) props.push(wrap({ ...values, [key]: option }));
        }
      } else if (p.type === 'boolean') {
        props.push(wrap({ ...values, [key]: !values[key] }));
      } else if (p.type === 'group') {
        vary(p.params, values[key], (v) => wrap({ ...values, [key]: v }));
      }
    }
  };
  for (const group of CanvasAxisGroups) {
    const param = (Canvas3DParams as any)[group] as PD.Group<any>;
    vary(param.params, defaults[group], (v) => ({ [group]: v }));
  }
  vary(
    { cameraFog: Canvas3DParams.cameraFog, transparentBackground: Canvas3DParams.transparentBackground },
    defaults,
    (v) => ({ cameraFog: v.cameraFog, transparentBackground: v.transparentBackground }),
  );
  // what a number turns on, and options that select programs together
  props.push(
    { renderer: { ...defaults.renderer, dimStrength: 0.5 } },
    {
      illumination: { ...defaults.illumination, enabled: true },
      multiSample: { ...defaults.multiSample, mode: 'on' },
    },
    // sharpening after FXAA and alone (SMAA, the default, needs images: not in Node)
    ...(['fxaa', 'off'] as const).map((name) => ({
      postprocessing: {
        ...defaults.postprocessing,
        antialiasing: { name, params: mapped('antialiasing').map(name).defaultValue },
        sharpening: { name: 'on', params: mapped('sharpening').map('on').defaultValue },
      },
    })),
  );
  const objects = unique(sources.flatMap((s) => s.list.objects));
  return { name: 'canvas axes', list: { objects }, props: props.flatMap(withTransparencyModes) };
}

/** The sources, each a case of its own (`compileSource`), in the order they are built */
export const SourceNames = ['presets', 'mvs', 'mesoscale', 'define axes', 'canvas axes'] as const;
export type SourceName = (typeof SourceNames)[number];

/**
 * The sources on `plugin`: the presets (built: their structures stay in the state, the real scene of `prepareFrames`),
 * the MolViewSpec documents, the mesoscale scenes, then the axes from them. Needs the MolViewSpec runtime in the
 * plugin, and a 2D canvas for the labels' glyphs
 */
export async function buildSources(plugin: PluginContext, read: ReadFile): Promise<VariantSource[]> {
  const presets = await getPresetSource(plugin, read);
  const mesoscale = getMesoscaleSource(plugin, getStructureRef(plugin));
  const mvs = await getMVSSource(plugin, read);
  const sources: VariantSource[] = [presets, mvs, mesoscale];
  sources.push(getDefineAxisSource(sources), getCanvasAxisSource(sources));
  await plugin.canvas3d!.graphics.resources.pending.whenReady();
  return sources;
}

/** The protein the presets loaded, by plugin: the parent of the mesoscale spacefills and the real scene's additions */
const Structures = new WeakMap<PluginContext, string>();
const getStructureRef = (plugin: PluginContext) => {
  const ref = Structures.get(plugin);
  if (!ref) throw new Error('no structure: the presets source is built first');
  return ref;
};

/** Plugins with the representations of `prepareFrames` added */
const withPassScene = new WeakSet<PluginContext>();

/**
 * Prepares the frames of the real scene (the presets' representations, and a transparent emissive one and a GPU surface
 * with smoothed colors, added once) with each of `props` set on the canvas, in its transparency mode: precompile plans
 * passes with stand-ins of the scene's groups, and without Hi-Z. Twice: what loads in the background (SMAA's lookup
 * images) is planned once loaded; the second time also drawing a frame on WebGPU, which compiles what preparing missed
 * (WebGL links every program preparing asks for). Returns what `count` was after the first and the second time
 */
export async function prepareFrames(
  plugin: PluginContext,
  props: readonly Canvas3DPrecompileProps[],
  count: () => number,
) {
  const canvas = plugin.canvas3d!;
  const { pending } = canvas.graphics.resources;
  if (!withPassScene.has(plugin)) {
    withPassScene.add(plugin);
    const structure = getStructureRef(plugin);
    await plugin.builders.structure.representation.addRepresentation(structure, {
      type: 'cartoon',
      typeParams: { alpha: 0.5, emissive: 0.5 },
    });
    // colors smoothed on the GPU (WebGL's are programs, run as the surface is built; WebGPU's kernels prepared apart)
    await plugin.builders.structure.representation.addRepresentation(structure, {
      type: 'gaussian-surface',
      typeParams: {
        tryUseGpu: true,
        // at most 1: WebGL computes coarser surfaces on the CPU
        quality: 'custom',
        resolution: 0.5,
        smoothColors: { name: 'on', params: { resolutionFactor: 2, sampleStride: 3 } },
      },
      color: 'chain-id',
    });
  }
  const byMode = new Map<string, Canvas3DPrecompileProps[]>();
  for (const { transparency = 'blended', ...p } of props) {
    const list = byMode.get(transparency) ?? [];
    if (!list.some((q) => JSON.stringify(q) === JSON.stringify(p))) list.push(p);
    byMode.set(transparency, list);
  }
  const initial = canvas.props;
  canvas.commit(true);
  // the viewport synced to the passes (a headless canvas is never resized): Hi-Z builds its levels for it
  canvas.handleResize();
  canvas.requestCameraReset();
  canvas.tick(now(), { isSynchronous: true });
  const draw = isWebGPUContext(canvas.graphics);
  let prepared = 0;
  for (let round = 0; round < 2; ++round) {
    for (const [transparency, list] of byMode) {
      plugin.canvas3dContext!.setProps({ transparency: transparency as TransparencyMode });
      canvas.commit(true);
      for (const p of list) {
        canvas.setProps({ ...initial, ...p });
        await canvas.prepare({ interaction: true });
        if (round === 1 && draw) {
          await pending.whenReady();
          canvas.requestDraw();
          canvas.tick(now(), { isSynchronous: true });
        }
      }
      await pending.whenReady();
    }
    if (round === 0) prepared = count();
  }
  canvas.setProps(initial);
  return { prepared, drawn: count(), props: [...byMode.values()].reduce((n, l) => n + l.length, 0) };
}

/** Compiles one source: its precompile list with its props, then its props as frames of the real scene */
export async function compileSource(plugin: PluginContext, source: VariantSource, count: () => number) {
  const before = count();
  await plugin.canvas3d!.precompile(source.list, source.props);
  await plugin.canvas3d!.graphics.resources.pending.whenReady();
  const precompiled = count();
  const frames = await prepareFrames(plugin, source.props, count);
  const n = source.list.objects.length;
  return (
    `${source.name}: ${n} signatures, ${source.props.length} props; precompiled ${precompiled - before}, ` +
    `prepared ${frames.prepared - precompiled}, drawn ${frames.drawn - frames.prepared}`
  );
}

/** The progress of a walk (`walkVariants`): what was compiled (the caller's count) by the end of each phase */
export interface VariantWalk {
  readonly sources: readonly VariantSource[];
  readonly summary: readonly string[];
  /** building the sources (the presets' frames) */
  readonly built: number;
  /** `Canvas3D.precompile` of each source */
  readonly precompiled: number;
  /** `prepare` of the real scene under every props */
  readonly prepared: number;
  /** a second `prepare` and a frame drawn under every props: what preparing missed */
  readonly drawn: number;
}

/**
 * Compiles every source (for the coverage of every program): builds them, precompiles each, then prepares the frames
 * of every props (`prepareFrames`). `count` is the caller's record of what was compiled so far (pipelines, programs)
 */
export async function walkVariants(plugin: PluginContext, read: ReadFile, count: () => number): Promise<VariantWalk> {
  const canvas = plugin.canvas3d!;
  const sources = await buildSources(plugin, read);
  const built = count();
  const summary: string[] = [];
  for (const source of sources) {
    const before = count();
    await canvas.precompile(source.list, source.props);
    await canvas.graphics.resources.pending.whenReady();
    const n = source.list.objects.length;
    summary.push(`${source.name}: ${n} signatures, ${source.props.length} props, ${count() - before} compiled`);
  }
  const precompiled = count();
  const frames = await prepareFrames(
    plugin,
    sources.flatMap((s) => s.props),
    count,
  );
  summary.push(
    `prepared frames: ${frames.props} props, ${frames.prepared - precompiled} compiled; ` +
      `drawn: ${frames.drawn - frames.prepared} more`,
  );
  return { sources, summary, built, precompiled, prepared: frames.prepared, drawn: frames.drawn };
}

/** What a graphics context compiled while recorded (`recordCompiles`) */
export interface CompileRecord {
  /** what was compiled so far (WebGPU: pipelines and kernel pipelines; WebGL: programs asked for), for `walkVariants` */
  readonly count: () => number;
  /** the registered programs (and kernels) that were compiled, among the first `count` compiled from `from` */
  readonly compiledIn: (from: number, to: number) => Set<string>;
  /** the registered programs (and kernels) never compiled */
  readonly uncovered: () => string[];
  /**
   * The compilers' errors and warnings, by message with the programs they were in: WebGPU, Dawn's (Tint's) messages
   * of every shader module; WebGL, the compile and link logs of failures
   */
  readonly problems: () => Promise<string[]>;
  /** puts the context's functions back */
  readonly stop: () => void;
}

/** Records what `graphics` compiles from now on, see `CompileRecord` */
export function recordCompiles(graphics: GraphicsContext): CompileRecord {
  return isWebGPUContext(graphics) ? recordWebGPU(graphics) : recordWebGL(graphics as WebGLContext);
}

/** problems by message, with the programs they were in */
function createProblems() {
  const problems = new Map<string, Set<string>>();
  return {
    add: (message: string, program: string) => {
      const programs = problems.get(message) ?? new Set<string>();
      programs.add(program);
      problems.set(message, programs);
    },
    list: () => [...problems].map(([message, programs]) => `${message} (${[...programs].sort().join(', ')})`),
  };
}

function recordWebGPU(graphics: WebGPUContext): CompileRecord {
  const { pending } = graphics.resources;
  const { device } = graphics;
  const problems = createProblems();
  const compiled: Promise<void>[] = [];
  const createShaderModule = device.createShaderModule;
  device.createShaderModule = (descriptor) => {
    const module = createShaderModule.call(device, descriptor);
    compiled.push(
      module.getCompilationInfo().then((info) => {
        const lines = descriptor.code.split('\n');
        for (const m of info.messages) {
          if (m.type === 'info') continue;
          const line = lines[m.lineNum - 1]?.trim() ?? '';
          problems.add(`${m.type}: ${m.message.split('\n')[0]} :: ${line}`, descriptor.label ?? '');
        }
      }),
    );
    return module;
  };
  // every label (the tracker keeps the recent ones only): `pipeline <program> <variant> <key>` (names may hold
  // spaces), `kernel pipeline <kernel>`
  const labels: string[] = [];
  const original = pending.track;
  const track = original.bind(pending);
  pending.track = (promise, label) => {
    if (label.startsWith('pipeline ') || label.startsWith('kernel pipeline ')) labels.push(label);
    return track(promise, label);
  };
  const names = [...getWgslProgramNames(), ...getWgslKernelNames()];
  const nameOf = (label: string) =>
    label.startsWith('kernel pipeline ')
      ? label.slice(16)
      : getWgslProgramNames().find((name) => label.startsWith(`pipeline ${name} `));
  const compiledIn = (from: number, to: number) => {
    const set = new Set<string>();
    for (const label of labels.slice(from, to)) {
      const name = nameOf(label);
      if (name) set.add(name);
    }
    return set;
  };
  return {
    count: () => labels.length,
    compiledIn,
    uncovered: () => {
      const compiled = compiledIn(0, labels.length);
      return names.filter((name) => !compiled.has(name));
    },
    problems: async () => {
      await Promise.all(compiled);
      return problems.list();
    },
    stop: () => {
      device.createShaderModule = createShaderModule;
      pending.track = original;
    },
  };
}

function recordWebGL(graphics: WebGLContext): CompileRecord {
  const { resources, gl } = graphics;
  const problems = createProblems();
  // every program asked for, by name (a cached one too)
  const programs: string[] = [];
  let current = '';
  const program = resources.program;
  resources.program = (defineValues, shaderCode, schema) => {
    programs.push(shaderCode.name);
    current = shaderCode.name;
    return program(defineValues, shaderCode, schema);
  };
  // the status of a shader or program is queried later, not right after compiling or linking: a query waits for the
  // driver, which would compile and link one after the other instead of in parallel (`KHR_parallel_shader_compile`);
  // before it is deleted, or when the problems are read
  const shaders = new Map<WebGLShader, string>();
  const linked = new Map<WebGLProgram, string>();
  const checkShader = (shader: WebGLShader, name: string) => {
    shaders.delete(shader);
    if (gl.getShaderParameter(shader, gl.COMPILE_STATUS)) return;
    const text = gl.getShaderSource(shader) ?? '';
    // the first error of the log, with the line it points at (`ERROR: 0:<line>: ...`), and the variant's defines
    const log = gl.getShaderInfoLog(shader)?.trim().split('\n')[0] ?? '';
    const line = Number(/^ERROR: \d+:(\d+):/.exec(log)?.[1]);
    const lines = text.split('\n');
    const source = lines[line - 1]?.trim() ?? '';
    const defines = lines
      .filter((l) => l.startsWith('#define d'))
      .map((l) => l.slice(8))
      .join(' ');
    problems.add(`compile: ${log} :: ${source} [${defines}]`, name);
  };
  const checkProgram = (p: WebGLProgram, name: string) => {
    linked.delete(p);
    if (!gl.getProgramParameter(p, gl.LINK_STATUS)) problems.add(`link: ${gl.getProgramInfoLog(p)?.trim()}`, name);
  };
  const compileShader = gl.compileShader;
  gl.compileShader = (shader) => {
    compileShader.call(gl, shader);
    shaders.set(shader, current);
  };
  const linkProgram = gl.linkProgram;
  gl.linkProgram = (p) => {
    linkProgram.call(gl, p);
    linked.set(p, current);
  };
  const deleteShader = gl.deleteShader;
  gl.deleteShader = (shader) => {
    const name = shader ? shaders.get(shader) : undefined;
    if (shader && name !== undefined) checkShader(shader, name);
    deleteShader.call(gl, shader);
  };
  const deleteProgram = gl.deleteProgram;
  gl.deleteProgram = (p) => {
    const name = p ? linked.get(p) : undefined;
    if (p && name !== undefined) checkProgram(p, name);
    deleteProgram.call(gl, p);
  };
  const compiledIn = (from: number, to: number) => new Set(programs.slice(from, to));
  return {
    count: () => programs.length,
    compiledIn,
    uncovered: () => {
      const compiled = compiledIn(0, programs.length);
      return getGlslProgramNames().filter((name) => !compiled.has(name));
    },
    problems: async () => {
      for (const [shader, name] of [...shaders]) checkShader(shader, name);
      for (const [p, name] of [...linked]) checkProgram(p, name);
      return problems.list();
    },
    stop: () => {
      resources.program = program;
      gl.compileShader = compileShader;
      gl.linkProgram = linkProgram;
      gl.deleteShader = deleteShader;
      gl.deleteProgram = deleteProgram;
    },
  };
}

/**
 * Prepares every compute kernel of the WebGPU backend the device can run, each request list in each variant (they
 * compile as they first run otherwise); nothing on WebGL, whose GPU paths are programs
 */
export async function prepareAllKernels(graphics: GraphicsContext) {
  if (!isWebGPUContext(graphics) || !graphics.compute) return;
  const { compute } = graphics;
  const requests = [
    ...GaussianDensityKernels,
    ...getColorSmoothingKernels('group'),
    ...getColorSmoothingKernels('groupInstance'),
    ...getIsosurfaceKernels({ kind: 'buffer' }),
    ...(['red', 'alpha'] as const).flatMap((channel) =>
      [false, true].flatMap((packedGroup) => getIsosurfaceKernels({ kind: 'texture', channel, packedGroup })),
    ),
    ...PrefixSumKernels,
    ...MergedCullKernels,
  ];
  // those the device can run (the core limits allow fewer storage buffers a stage), as the GPU paths choose
  await compute.prepareKernels(requests.filter((r) => compute.isKernelSupported(r.code, r.schema)));
  await graphics.resources.pending.whenReady();
}
