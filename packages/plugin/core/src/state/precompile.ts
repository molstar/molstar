/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { type StateTransform, type StateTransformer, StateTree } from '@molstar/core/state';
import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import type { Canvas3DPrecompileProps } from '@molstar/graphics/canvas3d/canvas3d';
import type { GraphicsContext } from '@molstar/graphics/gl/context';
import type { ObjectSignature, PrecompileList } from '@molstar/graphics/gl/precompile';
import type { GeometryUtils } from '@molstar/graphics/geo/geometry/geometry';
import {
  getColorThemeSignature,
  getColorTypes,
  getGeometrySignatures,
  getSizeThemeSignature,
  NoSignatureLayers,
  type SignatureContext,
  uniqueSignatures,
} from '@molstar/graphics/repr/signatures';
import type {
  Representation,
  RepresentationContext,
  RepresentationProvider,
} from '@molstar/graphics/repr/representation';
import { ColorTheme } from '@molstar/graphics/theme/color';
import { SizeTheme } from '@molstar/graphics/theme/size';
import type { ThemeRegistryContext } from '@molstar/graphics/theme/theme';
import type { PluginContext } from '../context.js';
import type { PluginState } from '../state.js';

/** What a transform's signatures are derived in: the plugin (registries), the context and the transform's place */
export interface TransformSignatureContext {
  readonly plugin: PluginContext;
  readonly graphics?: GraphicsContext;
  /** the layers set on the transform's representation by its children (overpaint, transparency, ...) */
  readonly layers: SignatureContext['layers'];
  readonly tree: StateTree;
  readonly transform: StateTransform;
}

/** The signatures of the objects a transform with `params` may create, without its data, see `getGeometrySignatures` */
export type TransformSignatures = (params: any, ctx: TransformSignatureContext) => ObjectSignature[];

type Layer = keyof SignatureContext['layers'];

const Signatures = new Map<string, TransformSignatures>();
const Layers = new Map<string, Layer>();

/** Registers how the signatures of `transformer`'s objects are derived from its params */
export function registerTransformSignatures(transformer: StateTransformer, signatures: TransformSignatures) {
  Signatures.set(transformer.id, signatures);
}

/** Registers `transformer` as one that sets `layer` on its parent representation */
export function registerTransformLayer(transformer: StateTransformer, layer: Layer) {
  Layers.set(transformer.id, layer);
}

function getLayers(tree: StateTree, ref: StateTransform.Ref): SignatureContext['layers'] {
  const children = tree.children.get(ref);
  if (!children || children.size === 0) return NoSignatureLayers;
  const layers = { ...NoSignatureLayers };
  children.forEach((c) => {
    const layer = Layers.get(tree.transforms.get(c!)!.transformer.id);
    if (layer) layers[layer] = true;
  });
  return layers;
}

function getParams(transform: StateTransform, plugin: PluginContext) {
  if (transform.params) return transform.params;
  const { params } = transform.transformer.definition;
  return params ? PD.getDefaultValues(params(undefined, plugin)) : {};
}

/**
 * The signatures of the objects the representations of `tree` may create, from their params alone (no data, nothing
 * built): what `Canvas3D.precompile` prepares before the tree is applied. Transforms without registered signatures
 * (`registerTransformSignatures`) are left out; their objects compile when first drawn
 */
export function getTreePrecompileList(
  plugin: PluginContext,
  tree: StateTree,
  cache?: Map<string, ObjectSignature[]>,
): PrecompileList {
  const graphics = plugin.canvas3d?.graphics;
  const signatures: ObjectSignature[] = [];
  StateTree.doPreOrder(tree, tree.root, undefined, (transform) => {
    const derive = Signatures.get(transform.transformer.id);
    if (!derive) return;
    try {
      const layers = getLayers(tree, transform.ref);
      const params = getParams(transform, plugin);
      // what the signatures follow from (shapes also from their provider, the parent)
      const parent = tree.transforms.get(transform.parent);
      const key =
        cache && JSON.stringify([transform.transformer.id, params, layers, parent?.transformer.id, parent?.params]);
      let derived = key !== undefined ? cache!.get(key) : undefined;
      if (!derived) {
        derived = derive(params, { plugin, graphics, layers, tree, transform });
        if (key !== undefined) cache!.set(key, derived);
      }
      signatures.push(...derived);
    } catch (e) {
      console.warn(`Deriving the signatures of '${transform.transformer.id}' failed`, e);
    }
  });
  return { objects: uniqueSignatures(signatures) };
}

/**
 * The list of `snapshots` (`getTreePrecompileList` of each state tree, and of the frames of their transitions, e.g. an
 * opacity animated) and the props their frames draw with, for `Canvas3D.precompile`
 */
export function getSnapshotsPrecompileList(plugin: PluginContext, snapshots: readonly PluginState.Snapshot[]) {
  const signatures: ObjectSignature[] = [];
  // snapshots and their frames often share transforms and props: each distinct ones once
  const cache = new Map<string, ObjectSignature[]>();
  const props = new Map<string, Canvas3DPrecompileProps>();
  const add = (
    data: { tree: StateTree.Serialized } | undefined,
    canvasProps: object | undefined,
    transparency?: string,
  ) => {
    if (data) signatures.push(...getTreePrecompileList(plugin, StateTree.fromJSON(data.tree), cache).objects);
    const p = { ...canvasProps, ...(transparency ? { transparency } : {}) } as Canvas3DPrecompileProps;
    props.set(JSON.stringify(p), p);
  };
  for (const snapshot of snapshots) {
    const transparency = snapshot.canvas3dContext?.props?.transparency;
    add(snapshot.data, snapshot.canvas3d?.props, transparency);
    for (const frame of snapshot.transition?.frames ?? []) {
      add(frame.data, frame.canvas3d?.props ?? snapshot.canvas3d?.props, transparency);
    }
  }
  return { list: { objects: uniqueSignatures(signatures) } as PrecompileList, props: [...props.values()] };
}

/**
 * Prepares what playing `snapshots` draws (`getSnapshotsPrecompileList`, `Canvas3D.precompile`), e.g. when a state or a
 * story is loaded; meant to run in the background (pending work), resolves when ready
 */
export function precompileSnapshots(plugin: PluginContext, snapshots: readonly PluginState.Snapshot[]): Promise<void> {
  const canvas = plugin.canvas3d;
  if (!canvas || snapshots.length === 0) return Promise.resolve();
  try {
    const { list, props } = getSnapshotsPrecompileList(plugin, snapshots);
    return canvas.precompile(list, props);
  } catch (e) {
    console.warn('Preparing the snapshots failed', e);
    return Promise.resolve();
  }
}

const AnyColorTypes: SignatureContext['color']['types'] = ['uniform', 'instance', 'group', 'groupInstance'];

function getColorSignature(themes: ThemeRegistryContext, theme: { name: string; params: any } | undefined) {
  const provider = theme && themes.colorThemeRegistry.get(theme.name);
  if (!provider) return getColorThemeSignature(ColorTheme.Empty);
  try {
    // without data, as the theme legends are made: themes that need data for their granularity name it
    const props = { ...provider.defaultValues, ...theme.params };
    const signature = getColorThemeSignature(provider.factory({}, props));
    const granularities = provider.getSignatureGranularities?.(props);
    return granularities ? { ...signature, types: granularities.flatMap(getColorTypes) } : signature;
  } catch {
    return { types: AnyColorTypes, palette: false, preferSmoothing: true };
  }
}

function getSizeSignature(themes: ThemeRegistryContext, theme: { name: string; params: any } | undefined) {
  const provider = theme && themes.sizeThemeRegistry.get(theme.name);
  if (!provider) return getSizeThemeSignature(SizeTheme.Empty);
  try {
    return getSizeThemeSignature(provider.factory({}, { ...provider.defaultValues, ...theme.params }));
  } catch {
    return { types: ['uniform', 'instance', 'group', 'groupInstance'] as const };
  }
}

/**
 * The signatures of a representation from a registry with `params` (type, color and size themes, as the representation
 * transforms take them), created without data
 */
export function getRepresentationSignatures(
  provider: RepresentationProvider<any, any, any> | undefined,
  themes: ThemeRegistryContext,
  params: { type: { name: string; params?: any }; colorTheme?: { name: string; params: any }; sizeTheme?: any },
  ctx: TransformSignatureContext,
): ObjectSignature[] {
  if (!provider) return [];
  const repr = provider.factory({ graphics: ctx.graphics, ...themes }, provider.getParams);
  try {
    // an undefined param is the default, as in a state read from JSON (MVS writes `clip: undefined` without clip nodes)
    const props: { [k: string]: unknown } = { ...provider.defaultValues };
    for (const [k, v] of Object.entries(params.type.params ?? {})) if (v !== undefined) props[k] = v;
    return (
      repr.getSignatures?.(props, {
        graphics: ctx.graphics,
        color: getColorSignature(themes, params.colorTheme),
        size: getSizeSignature(themes, params.sizeTheme),
        layers: ctx.layers,
      }) ?? []
    );
  } finally {
    repr.destroy();
  }
}

/** The geometries of the shapes a transform provides (`SO.Shape.Provider`) with `params` */
export type ShapeGeometries = (params: any) => GeometryUtils<any>[];

const ShapeGeometry = new Map<string, ShapeGeometries>();

/** Registers the geometries of the shapes `transformer` provides, for the signatures of their representations */
export function registerShapeGeometry(transformer: StateTransformer, geometries: ShapeGeometries) {
  ShapeGeometry.set(transformer.id, geometries);
}

/**
 * The signatures of a shape representation with `params` of the shapes `provider` (its parent transform) gives, colored
 * and sized by shape group (see `ShapeRepresentation`)
 */
export function getShapeSignatures(
  provider: StateTransform | undefined,
  params: any,
  ctx: TransformSignatureContext,
): ObjectSignature[] {
  const geometries = provider && ShapeGeometry.get(provider.transformer.id);
  if (!geometries) return [];
  return geometries(getParams(provider, ctx.plugin)).flatMap((utils) =>
    getGeometrySignatures(utils, { ...PD.getDefaultValues(utils.Params), ...params }, getShapeContext(ctx)),
  );
}

/** Shapes are colored and sized by shape group (see `ShapeRepresentation`) */
function getShapeContext(ctx: TransformSignatureContext): SignatureContext {
  return {
    graphics: ctx.graphics,
    color: { types: ['groupInstance'], palette: false, preferSmoothing: false },
    size: { types: ['groupInstance'] },
    layers: ctx.layers,
  };
}

/**
 * The signatures of a representation of shapes (e.g. a measurement's) that `create` makes, with `params`, created
 * without data
 */
export function getShapeRepresentationSignatures(
  create: (ctx: RepresentationContext) => Representation<any, any>,
  params: any,
  ctx: TransformSignatureContext,
): ObjectSignature[] {
  const repr = create({ graphics: ctx.graphics, ...ctx.plugin.representation.structure.themes });
  try {
    return repr.getSignatures?.({ ...repr.props, ...params }, getShapeContext(ctx)) ?? [];
  } finally {
    repr.destroy();
  }
}
