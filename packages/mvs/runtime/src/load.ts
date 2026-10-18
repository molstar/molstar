/**
 * Copyright (c) 2023-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Adam Midlik <midlik@gmail.com>
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Aliaksei Chareshneu <chareshneu.tech@gmail.com>
 * @author Marek Eibel <108737044+kerrambit@users.noreply.github.com>
 */

import { PluginStateSnapshotManager } from '@molstar/plugin/state/manager/snapshots';
import { PluginStateObject } from '@molstar/plugin/state/objects';
import { Download } from '@molstar/plugin/state/transforms/data/fetch';
import { ParseCcp4 } from '@molstar/plugin/state/formats/volume/ccp4';
import { ParseCif } from '@molstar/plugin/state/formats/cif';
import { ParseDx } from '@molstar/plugin/state/formats/volume/dx';
import { ParseObj } from '@molstar/plugin/state/formats/shape/obj';
import { ParsePly } from '@molstar/plugin/state/formats/shape/ply';
import { ParsePrmtop } from '@molstar/plugin/state/formats/topology/prmtop';
import { ParsePsf } from '@molstar/plugin/state/formats/topology/psf';
import { ParseTop } from '@molstar/plugin/state/formats/topology/top';
import { ParseVtp } from '@molstar/plugin/state/formats/shape/vtp';
import { CoordinatesFromDcd } from '@molstar/plugin/state/formats/coordinates/dcd';
import { CoordinatesFromLammpstraj } from '@molstar/plugin/state/formats/coordinates/lammps';
import { CoordinatesFromNctraj } from '@molstar/plugin/state/formats/coordinates/nctraj';
import { CoordinatesFromTrr } from '@molstar/plugin/state/formats/coordinates/trr';
import { CoordinatesFromXtc } from '@molstar/plugin/state/formats/coordinates/xtc';
import {
  CustomModelProperties,
  CustomStructureProperties,
  ModelFromTrajectory,
  StructureFromModel,
} from '@molstar/plugin/state/transforms/structure/hierarchy';
import type { StructureComponentParams } from '@molstar/plugin/state/helpers/structure-component';
import { StructureComponent } from '@molstar/plugin/state/transforms/structure/selection';
import { TopologyFromPrmtop } from '@molstar/plugin/state/formats/topology/prmtop';
import { TopologyFromPsf } from '@molstar/plugin/state/formats/topology/psf';
import { TopologyFromTop } from '@molstar/plugin/state/formats/topology/top';
import { TrajectoryFromGRO } from '@molstar/plugin/state/formats/trajectory/gro';
import { TrajectoryFromLammpsTrajData } from '@molstar/plugin/state/formats/trajectory/lammps';
import { TrajectoryFromMmCif } from '@molstar/plugin/state/formats/trajectory/mmcif';
import { TrajectoryFromMOL } from '@molstar/plugin/state/formats/trajectory/mol';
import { TrajectoryFromMOL2 } from '@molstar/plugin/state/formats/trajectory/mol2';
import { TrajectoryFromPDB } from '@molstar/plugin/state/formats/trajectory/pdb';
import { TrajectoryFromSDF } from '@molstar/plugin/state/formats/trajectory/sdf';
import { TrajectoryFromXYZ } from '@molstar/plugin/state/formats/trajectory/xyz';
import { ShapeRepresentation3D } from '@molstar/plugin/state/transforms/shape/representation';
import { StructureRepresentation3D } from '@molstar/plugin/state/transforms/structure/representation';
import { VolumeRepresentation3D } from '@molstar/plugin/state/transforms/volume/representation';
import { ShapeFromObj } from '@molstar/plugin/state/formats/shape/obj';
import { ShapeFromPly } from '@molstar/plugin/state/formats/shape/ply';
import { ShapeFromVtp } from '@molstar/plugin/state/formats/shape/vtp';
import { VolumeFromCcp4 } from '@molstar/plugin/state/formats/volume/ccp4';
import { VolumeFromDensityServerCif } from '@molstar/plugin/state/formats/volume/density-server';
import { VolumeFromDx } from '@molstar/plugin/state/formats/volume/dx';
import { PluginCommands } from '@molstar/plugin/commands';
import { PluginContext } from '@molstar/plugin/context';
import { PluginState } from '@molstar/plugin/state';
import { type State, StateObjectSelector, StateTree } from '@molstar/core/state';
import { RuntimeContext, Task } from '@molstar/core/task';
import { Clip } from '@molstar/core/util/clip';
import { ScratchMemory } from '@molstar/core/util/int-buffer';
import { SchemaSelectionScratch } from '@molstar/model/model/structure/structure/element/schema-selection/scratch';
import { MVS_BEHAVIOR_ID } from './behavior-id.js';
import { createPluginStateSnapshotCamera, modifyCanvasProps, resetCanvasProps } from '@molstar/mvs/camera';
import { MVSAnnotationsProvider } from '@molstar/mvs/components/annotation-prop';
import { MVSAnnotationStructureComponent } from '@molstar/mvs/components/annotation-structure-component';
import { MVSAnnotationTooltipsProvider } from '@molstar/mvs/components/annotation-tooltips-prop';
import {
  type CustomLabelProps,
  CustomLabelRepresentationProvider,
} from '@molstar/mvs/components/custom-label/representation';
import { CustomTooltipsProvider } from '@molstar/mvs/components/custom-tooltips-prop';
import { type IsMVSModelProps, IsMVSModelProvider } from '@molstar/mvs/components/is-mvs-model-prop';
import {
  getPrimitiveStructureRefs,
  MVSBuildPrimitiveShape,
  MVSDownloadPrimitiveData,
  MVSInlinePrimitiveData,
  MVSShapeRepresentation3D,
} from '@molstar/mvs/components/primitives';
import { MVSTrajectoryWithCoordinates } from '@molstar/mvs/components/trajectory';
import { generateStateTransition } from '@molstar/mvs/helpers/animation';
import { type ShapeFormat, shapeRepresentationProps, shapeTransforms } from '@molstar/mvs/helpers/load-shape';
import { IsHiddenCustomStateExtension } from '@molstar/mvs/load-extensions/is-hidden-custom-state';
import { NonCovalentInteractionsExtension } from '@molstar/mvs/load-extensions/non-covalent-interactions';
import { VolumeStreamingExtension } from '@molstar/mvs/load-extensions/volume-streaming';
import {
  type LoadingActions,
  type LoadingExtension,
  loadTreeVirtual,
  SubtreeMemo,
  UpdateTarget,
} from '@molstar/mvs/load-generic';
import {
  type AnnotationFromSourceKind,
  type AnnotationFromUriKind,
  clippingForNode,
  collectAnnotationReferences,
  collectAnnotationTooltips,
  collectInlineLabels,
  collectInlineTooltips,
  colorThemeForNode,
  componentFromXProps,
  componentPropsFromSelector,
  isPhantomComponent,
  labelFromXProps,
  makeNearestReprMap,
  prettyNameFromSelector,
  representationProps,
  structureProps,
  transformAndInstantiateStructure,
  transformAndInstantiateVolume,
  volumeColorThemeForNode,
  volumeRepresentationProps,
} from '@molstar/mvs/load-helpers';
import { MVSData, type MVSData_States, type Snapshot, SnapshotMetadata } from '@molstar/mvs-builder/mvs-data';
import { type MVSAnimationNode, MVSAnimationSchema } from '@molstar/mvs-builder/tree/animation/animation-tree';
import { validateTree } from '@molstar/mvs-builder/tree/generic/tree-validation';
import {
  convertMvsToMolstar,
  convertMvsToMolstarMemo,
  createMolstarConversionMemo,
  mvsSanityCheck,
} from '@molstar/mvs/tree/molstar/conversion';
import {
  type MolstarNode,
  type MolstarNodeParams,
  type MolstarSubtree,
  type MolstarTree,
  MolstarTreeSchema,
} from '@molstar/mvs/tree/molstar/molstar-tree';
import { MVSTreeSchema } from '@molstar/mvs-builder/tree/mvs/mvs-tree';
import { validateMolQLTree } from '@molstar/mvs/runtime-validation';

export interface MVSLoadOptions {
  /** Add snapshots from MVS into current snapshot list, instead of replacing the list. */
  appendSnapshots?: boolean;
  /** Ignore any camera positioning from the MVS state and keep the current camera position instead, ignore any camera positioning when generating snapshots. */
  keepCamera?: boolean;
  /** Follow camera target position from the MVS state but keep the current camera direction and up. (`keepCamera` option overrides this) */
  keepCameraOrientation?: boolean;
  /** Specifies a set of MVS-loading extensions (not a part of standard MVS specification). If undefined, apply all builtin extensions. If `[]`, do not apply builtin extensions. */
  extensions?: MolstarLoadingExtension<any>[];
  /** Run some sanity checks and print potential issues to the console. */
  sanityChecks?: boolean;
  /** Base for resolving relative URLs/URIs. May itself be a relative URL (relative to the window URL). */
  sourceUrl?: string;
  /** The index of the snapshot to apply initially. Defaults to 0. */
  defaultSnapshotIndex?: number;

  doNotReportErrors?: boolean;
}

export function loadMVS(plugin: PluginContext, data: MVSData, options: MVSLoadOptions = {}) {
  const task = Task.create('Load MVS', (ctx) => _loadMVS(ctx, plugin, data, options));
  return plugin.runTask(task);
}

/** Load a MolViewSpec (MVS) state(s) into the Mol* plugin as plugin state snapshots. */
async function _loadMVS(ctx: RuntimeContext, plugin: PluginContext, data: MVSData, options: MVSLoadOptions = {}) {
  plugin.errorContext.clear('mvs');
  try {
    const mvsExtensionLoaded = plugin.state.hasBehavior(MVS_BEHAVIOR_ID);
    if (!mvsExtensionLoaded) throw new Error('MolViewSpec extension is not loaded.');

    // Stop any currently running audio
    plugin.managers.markdownExtensions.audio.dispose();

    // Reset canvas props to default so that modifyCanvasProps works as expected
    resetCanvasProps(plugin);

    // console.log(`MVS tree:\n${MVSData.toPrettyString(data)}`)
    const multiData: MVSData_States = data.kind === 'multiple' ? data : MVSData.stateToStates(data);
    const entries: PluginStateSnapshotManager.Entry[] = [];
    // (refs are the same in every snapshot and frame: their hashes are computed once per load, `load.md §2`)
    const loadOptions: MolstarTreeLoadOptions = { ...options, refHashCache: new Map() };
    for (let i = 0; i < multiData.snapshots.length; i++) {
      const snapshot = multiData.snapshots[i];
      const previousSnapshot = i > 0 ? multiData.snapshots[i - 1] : multiData.snapshots[multiData.snapshots.length - 1];
      validateTree(MVSTreeSchema, snapshot.root, 'MVS', plugin);
      validateMolQLTree(snapshot.root, 'MVS', plugin);
      if (snapshot.animation) {
        validateTree(MVSAnimationSchema, snapshot.animation, 'Animation', plugin);
        validateMolQLTree(snapshot.animation, 'Animation', plugin);
      }
      if (options.sanityChecks) mvsSanityCheck(snapshot.root);
      const molstarTree = convertMvsToMolstar(snapshot.root, options.sourceUrl);
      validateTree(MolstarTreeSchema, molstarTree, 'Converted Molstar', plugin);
      const entry = molstarTreeToEntry(
        plugin,
        molstarTree,
        snapshot.animation,
        { ...snapshot.metadata, previousTransitionDurationMs: previousSnapshot.metadata.transition_duration_ms },
        loadOptions,
      );
      await assignStateTransition(ctx, plugin, entry, snapshot, loadOptions, i, multiData.snapshots.length);
      entries.push({
        ...entry,
        _transientData: { sourceMvsSnapshot: snapshot },
      });

      if (ctx.shouldUpdate) {
        await ctx.update({ message: 'Loading MVS...', current: i, max: multiData.snapshots.length });
      }
    }
    if (!options.appendSnapshots) {
      plugin.managers.snapshot.clear();
    }
    for (const entry of entries) {
      plugin.managers.snapshot.add(entry);
    }

    if (entries.length > 0) {
      let indexToApply = options.defaultSnapshotIndex ?? 0;
      if (indexToApply < 0 || indexToApply >= entries.length) {
        indexToApply = 0;
      }
      await PluginCommands.State.Snapshots.Apply(plugin, { id: entries[indexToApply].snapshot.id });
    }
  } catch (err) {
    plugin.log.error(`${err}`);
    throw err;
  } finally {
    if (!options.doNotReportErrors) {
      for (const error of plugin.errorContext.get('mvs')) {
        plugin.log.warn(error);
        PluginCommands.Toast.Show(plugin, {
          title: 'Error',
          message: error,
          timeoutMs: 10000,
        });
      }
    }
    plugin.errorContext.clear('mvs');
    // (the buffers which the resolution of rows reuses between calls are not needed until the next load)
    ScratchMemory.release(SchemaSelectionScratch);
  }
}

async function assignStateTransition(
  ctx: RuntimeContext,
  plugin: PluginContext,
  parentEntry: PluginStateSnapshotManager.Entry,
  parent: Snapshot,
  options: MolstarTreeLoadOptions,
  snapshotIndex: number,
  snapshotCount: number,
) {
  const transitions = await generateStateTransition(ctx, parent, snapshotIndex, snapshotCount);
  if (!transitions?.frames.length) return;

  const animation: PluginState.StateTransition = {
    autoplay: !!transitions.tree.params?.autoplay,
    loop: !!transitions.tree.params?.loop,
    frames: [],
  };
  const includeCamera = !!transitions.tree.params?.include_camera;
  const includeCanvas = !!transitions.tree.params?.include_canvas;
  const parentData = parentEntry.snapshot.data!;
  // (a subtree that is the same object as in an earlier frame converts and loads as then, `load.md §4`; not with
  // camera or canvas, which the loading context collects from the whole tree)
  const memoize = !includeCamera && !includeCanvas;
  const conversionMemo = createMolstarConversionMemo(options.sourceUrl);
  const subtreeMemo = new SubtreeMemo<MolstarTree>(isReusableSubtree);

  let previousData: State.Snapshot | undefined;
  for (let i = 0; i < transitions.frames.length; i++) {
    const frame = transitions.frames[i];
    let molstarTree: MolstarTree;
    if (memoize) {
      const converted = convertMvsToMolstarMemo(frame[0], conversionMemo);
      molstarTree = converted.tree;
      subtreeMemo.units = converted.units;
    } else {
      molstarTree = convertMvsToMolstar(frame[0], options.sourceUrl);
    }
    // (a frame's transforms that equal the previous frame's or the snapshot's are those objects, `load.md §3`)
    const reuse = previousData ? [previousData.tree, parentData.tree] : [parentData.tree];
    let data: State.Snapshot;
    let camera: PluginState.Snapshot['camera'];
    let canvas3d: PluginState.Snapshot['canvas3d'];
    if (includeCamera || includeCanvas) {
      const entry = molstarTreeToEntry(
        plugin,
        molstarTree,
        parent.animation,
        { ...parent.metadata, previousTransitionDurationMs: transitions.frametimeMs },
        { ...options, reuse },
      );
      data = entry.snapshot.data!;
      camera = entry.snapshot.camera;
      canvas3d = entry.snapshot.canvas3d;
    } else {
      // (only the state tree: the camera and canvas of a frame are not kept without include_camera / include_canvas)
      data = loadTreeVirtual(plugin, molstarTree, MolstarLoadingActions, MolstarLoadingContext.create(), {
        replaceExisting: true,
        extensions: options.extensions ?? BuiltinLoadingExtensions,
        refHashCache: options.refHashCache,
        reuse,
        subtreeMemo: memoize ? subtreeMemo : undefined,
      }).data!;
    }
    previousData = data;

    animation.frames.push({
      durationInMs: frame[1],
      data,
      camera: includeCamera ? camera : undefined,
      canvas3d: includeCanvas ? canvas3d : undefined,
    });

    if (ctx.shouldUpdate) {
      await ctx.update({
        message: `Loading animation for snapshot ${snapshotIndex + 1}/${snapshotCount}...`,
        current: i + 1,
        max: transitions.frames.length,
      });
    }
  }

  parentEntry.snapshot.transition = animation;
}

/** Kinds whose loading reads outside their subtree (the nearest representation, `makeNearestReprMap`). */
const NonLocalKinds = new Set<string>(['label_from_uri', 'label_from_source']);

/** A unit subtree can be replayed unless loading it reads outside it (`load.md §4`). */
function isReusableSubtree(node: MolstarSubtree): boolean {
  if (NonLocalKinds.has(node.kind)) return false;
  for (const child of (node.children ?? []) as MolstarSubtree[]) {
    if (!isReusableSubtree(child)) return false;
  }
  return true;
}

interface MolstarTreeLoadOptions extends MVSLoadOptions {
  /** Ref hashes shared by the snapshots and frames of a load. */
  refHashCache?: Map<string, string>;
  /** Serialized trees whose equal transforms are reused (`loadTreeVirtual`). */
  reuse?: StateTree.Serialized[];
}

function molstarTreeToEntry(
  plugin: PluginContext,
  tree: MolstarTree,
  animation: MVSAnimationNode<'animation'> | undefined,
  metadata: SnapshotMetadata & { previousTransitionDurationMs?: number },
  options: MolstarTreeLoadOptions,
) {
  const context = MolstarLoadingContext.create();
  const snapshot = loadTreeVirtual(plugin, tree, MolstarLoadingActions, context, {
    replaceExisting: true,
    extensions: options?.extensions ?? BuiltinLoadingExtensions,
    refHashCache: options.refHashCache,
    reuse: options.reuse,
  });
  snapshot.canvas3d = {
    props: plugin.canvas3d ? modifyCanvasProps(plugin.canvas3d.props, context.canvas, animation) : undefined,
  };
  const snapshotDurationMs =
    metadata.duration_ms ?? metadata.linger_duration_ms ?? SnapshotMetadata.Defaults.duration_ms;
  const transitionParams = context.transition?.params;
  const transitionDurationMs = transitionParams?.duration_ms ?? metadata.previousTransitionDurationMs ?? 0;
  if (options?.keepCamera) {
    // do nothing
  } else if (options.keepCameraOrientation) {
    // load camera target, keep orientation
    snapshot.camera = createPluginStateSnapshotCamera(plugin, context, {
      incomingTransitionDurationMs: transitionDurationMs,
      ignoreCameraOrientation: true,
    });
  } else {
    // fully load camera
    snapshot.camera = createPluginStateSnapshotCamera(plugin, context, {
      incomingTransitionDurationMs: transitionDurationMs,
    });
  }
  if (snapshot.camera) {
    if (transitionParams?.easing) snapshot.camera.transitionEasing = transitionParams.easing;
    if (transitionParams?.trajectory) snapshot.camera.transitionTrajectory = transitionParams.trajectory;
  }
  snapshot.durationInMs = snapshotDurationMs + transitionDurationMs;
  snapshot.structureFocus = {}; // avoid structure focus persisting through states (causes weird behaviors, e.g. when turning on Volume Streaming)

  if (tree.custom?.molstar_on_load_markdown_commands) {
    snapshot.onLoadMarkdownCommands = tree.custom.molstar_on_load_markdown_commands;
  }

  const entryParams: PluginStateSnapshotManager.EntryParams = {
    key: metadata.key,
    name: metadata.title,
    description: metadata.description,
    descriptionFormat: metadata.description_format ?? 'markdown',
  };
  const entry: PluginStateSnapshotManager.Entry = PluginStateSnapshotManager.Entry(snapshot, entryParams);
  return entry;
}

/** Mutable context for loading a `MolstarTree`, available throughout the loading. */
export interface MolstarLoadingContext {
  /** Maps `*_from_[uri|source]` nodes to annotationId they should reference */
  annotationMap: Map<MolstarNode<AnnotationFromUriKind | AnnotationFromSourceKind>, string>;
  /** Maps each node (on 'structure' or lower level) to its nearest 'representation' node */
  nearestReprMap?: Map<MolstarNode, MolstarNode<'representation'>>;
  camera: {
    cameraParams?: MolstarNodeParams<'camera'>;
    focuses: { target: StateObjectSelector; params: MolstarNodeParams<'focus'> }[];
  };
  canvas?: MolstarNode<'canvas'>;
  transition?: MolstarNode<'transition'>;
}
export const MolstarLoadingContext = {
  create(): MolstarLoadingContext {
    return {
      annotationMap: new Map(),
      camera: { focuses: [] },
    };
  },
};

/** Loading actions for loading a `MolstarTree`, per node kind. */
const MolstarLoadingActions: LoadingActions<MolstarTree, MolstarLoadingContext> = {
  root(updateParent: UpdateTarget, node: MolstarNode<'root'>, context: MolstarLoadingContext): UpdateTarget {
    context.nearestReprMap = makeNearestReprMap(node);
    return updateParent;
  },
  download(updateParent: UpdateTarget, node: MolstarNode<'download'>): UpdateTarget {
    return UpdateTarget.apply(updateParent, Download, {
      url: node.params.url,
      isBinary: node.params.is_binary,
    });
  },
  parse(updateParent: UpdateTarget, node: MolstarNode<'parse'>): UpdateTarget | undefined {
    const format = node.params.format;
    switch (format) {
      case 'cif':
        return UpdateTarget.apply(updateParent, ParseCif, {});
      case 'pdb':
      case 'pdbqt':
      case 'gro':
      case 'xyz':
      case 'mol':
      case 'sdf':
      case 'mol2':
      case 'xtc':
      case 'lammpstrj':
      case 'dcd':
      case 'nctraj':
      case 'trr':
        return updateParent;
      case 'psf':
        return UpdateTarget.apply(updateParent, ParsePsf, {});
      case 'prmtop':
        return UpdateTarget.apply(updateParent, ParsePrmtop, {});
      case 'top':
        return UpdateTarget.apply(updateParent, ParseTop, {});
      case 'map':
        return UpdateTarget.apply(updateParent, ParseCcp4, {});
      case 'dx':
      case 'dxbin':
        return UpdateTarget.apply(updateParent, ParseDx, {});
      case 'vtp':
        return UpdateTarget.apply(updateParent, ParseVtp, {});
      case 'ply':
        return UpdateTarget.apply(updateParent, ParsePly, {});
      case 'obj':
        return UpdateTarget.apply(updateParent, ParseObj, {});
      default:
        console.error(`Unknown format in "parse" node: "${format}"`);
        return undefined;
    }
  },
  coordinates(updateParent: UpdateTarget, node: MolstarNode<'coordinates'>): UpdateTarget | undefined {
    const format = node.params.format;
    switch (format) {
      case 'nctraj':
        return UpdateTarget.apply(updateParent, CoordinatesFromNctraj);
      case 'dcd':
        return UpdateTarget.apply(updateParent, CoordinatesFromDcd);
      case 'trr':
        return UpdateTarget.apply(updateParent, CoordinatesFromTrr);
      case 'xtc':
        return UpdateTarget.apply(updateParent, CoordinatesFromXtc);
      case 'lammpstrj':
        return UpdateTarget.apply(updateParent, CoordinatesFromLammpstraj);
      default:
        console.error(`Unknown format in "coordinates" node: "${format}"`);
        return undefined;
    }
  },
  trajectory(updateParent: UpdateTarget, node: MolstarNode<'trajectory'>): UpdateTarget | undefined {
    const format = node.params.format;
    switch (format) {
      case 'cif':
        return UpdateTarget.apply(updateParent, TrajectoryFromMmCif, {
          blockHeader: node.params.block_header ?? '', // Must set to '' because just undefined would get overwritten by createDefaults
          blockIndex: node.params.block_index ?? undefined,
        });
      case 'pdb':
      case 'pdbqt':
        return UpdateTarget.apply(updateParent, TrajectoryFromPDB, { variant: format });
      case 'gro':
        return UpdateTarget.apply(updateParent, TrajectoryFromGRO);
      case 'xyz':
        return UpdateTarget.apply(updateParent, TrajectoryFromXYZ);
      case 'mol':
        return UpdateTarget.apply(updateParent, TrajectoryFromMOL);
      case 'sdf':
        return UpdateTarget.apply(updateParent, TrajectoryFromSDF);
      case 'mol2':
        return UpdateTarget.apply(updateParent, TrajectoryFromMOL2);
      case 'lammpstrj':
        return UpdateTarget.apply(updateParent, TrajectoryFromLammpsTrajData);
      default:
        console.error(`Unknown format in "trajectory" node: "${format}"`);
        return undefined;
    }
  },
  trajectory_with_coordinates(
    updateParent: UpdateTarget,
    node: MolstarNode<'trajectory_with_coordinates'>,
  ): UpdateTarget | undefined {
    const result = UpdateTarget.apply(updateParent, MVSTrajectoryWithCoordinates, {
      coordinatesRef: node.params.coordinates_ref,
    });
    return UpdateTarget.setMvsDependencies(result, [node.params.coordinates_ref]);
  },
  topology_with_coordinates(
    updateParent: UpdateTarget,
    node: MolstarNode<'topology_with_coordinates'>,
  ): UpdateTarget | undefined {
    let parsed: UpdateTarget;
    const format = node.params.format;
    switch (format) {
      case 'psf':
        parsed = UpdateTarget.apply(updateParent, TopologyFromPsf, {});
        break;
      case 'prmtop':
        parsed = UpdateTarget.apply(updateParent, TopologyFromPrmtop, {});
        break;
      case 'top':
        parsed = UpdateTarget.apply(updateParent, TopologyFromTop, {});
        break;
      default:
        console.error(`Unknown format in "topology_with_coordinates" node: "${format}"`);
        return undefined;
    }
    const result = UpdateTarget.apply(parsed, MVSTrajectoryWithCoordinates, {
      coordinatesRef: node.params.coordinates_ref,
    });
    return UpdateTarget.setMvsDependencies(result, [node.params.coordinates_ref]);
  },
  model(updateParent: UpdateTarget, node: MolstarSubtree<'model'>, context: MolstarLoadingContext): UpdateTarget {
    const annotations = collectAnnotationReferences(node, context);
    const model = UpdateTarget.apply(updateParent, ModelFromTrajectory, {
      modelIndex: node.params.model_index,
    });
    UpdateTarget.apply(model, CustomModelProperties, {
      properties: {
        [IsMVSModelProvider.descriptor.name]: { isMvs: true } satisfies IsMVSModelProps,
        [MVSAnnotationsProvider.descriptor.name]: { annotations },
      },
      // autoAttach not needed for these properties because they have isHidden:true (explicit autoAttach here would override default autoAttach of other properties!, e.g. sifts_sequence_mapping)
    });
    return model;
  },
  structure(
    updateParent: UpdateTarget,
    node: MolstarSubtree<'structure'>,
    context: MolstarLoadingContext,
  ): UpdateTarget {
    const props = structureProps(node);
    const struct = UpdateTarget.apply(updateParent, StructureFromModel, props);
    const transformed = transformAndInstantiateStructure(struct, node);
    const annotationTooltips = collectAnnotationTooltips(node, context);
    const inlineTooltips = collectInlineTooltips(node, context);
    UpdateTarget.apply(struct, CustomStructureProperties, {
      properties: {
        [MVSAnnotationTooltipsProvider.descriptor.name]: { tooltips: annotationTooltips },
        [CustomTooltipsProvider.descriptor.name]: { tooltips: inlineTooltips },
      },
      // autoAttach not needed for these properties because they have isHidden:true (explicit autoAttach here would override default autoAttach of other properties!)
    }); // CustomStructureProperties must be applied even when `annotationTooltips` and `inlineTooltips` are empty, otherwise tooltips would persists across MVS snapshots
    const inlineLabels = collectInlineLabels(node, context);
    if (inlineLabels.length > 0) {
      const nearestReprNode = context.nearestReprMap?.get(node);
      UpdateTarget.apply(struct, StructureRepresentation3D, {
        type: {
          name: CustomLabelRepresentationProvider.name,
          params: { items: inlineLabels } satisfies Partial<CustomLabelProps>,
        },
        colorTheme: colorThemeForNode(nearestReprNode, context),
      });
    }
    return transformed;
  },
  tooltip: undefined, // No action needed, already loaded in `structure`
  tooltip_from_uri: undefined, // No action needed, already loaded in `structure`
  tooltip_from_source: undefined, // No action needed, already loaded in `structure`
  component(updateParent: UpdateTarget, node: MolstarSubtree<'component'>): UpdateTarget | undefined {
    if (isPhantomComponent(node)) {
      return updateParent;
    }
    const selector = node.params.selector;
    const type = componentPropsFromSelector(selector);
    const label = prettyNameFromSelector(selector);
    // (rows are the expression `element-schema`, which resolves them through the indices of the model; a `Selector` has `annotation`, which a component does not)
    const component = UpdateTarget.apply(updateParent, StructureComponent, {
      type: type as StructureComponentParams['type'],
      label,
      nullIfEmpty: false,
    });
    return transformAndInstantiateStructure(component, node);
  },
  component_from_uri(
    updateParent: UpdateTarget,
    node: MolstarSubtree<'component_from_uri'>,
    context: MolstarLoadingContext,
  ): UpdateTarget | undefined {
    if (isPhantomComponent(node)) return undefined;
    const props = componentFromXProps(node, context);
    return transformAndInstantiateStructure(
      UpdateTarget.apply(updateParent, MVSAnnotationStructureComponent, props),
      node,
    );
  },
  component_from_source(
    updateParent: UpdateTarget,
    node: MolstarSubtree<'component_from_source'>,
    context: MolstarLoadingContext,
  ): UpdateTarget | undefined {
    if (isPhantomComponent(node)) return undefined;
    const props = componentFromXProps(node, context);
    return transformAndInstantiateStructure(
      UpdateTarget.apply(updateParent, MVSAnnotationStructureComponent, props),
      node,
    );
  },
  representation(
    updateParent: UpdateTarget,
    node: MolstarNode<'representation'>,
    context: MolstarLoadingContext,
  ): UpdateTarget {
    return UpdateTarget.apply(updateParent, StructureRepresentation3D, {
      ...representationProps(node),
      colorTheme: colorThemeForNode(node, context),
    });
  },
  volume(updateParent: UpdateTarget, node: MolstarNode<'volume'>): UpdateTarget | undefined {
    let volume: UpdateTarget;
    if (updateParent.transformer?.definition.to.includes(PluginStateObject.Format.Ccp4)) {
      volume = UpdateTarget.apply(updateParent, VolumeFromCcp4, {});
    } else if (updateParent.transformer?.definition.to.includes(PluginStateObject.Format.Dx)) {
      volume = UpdateTarget.apply(updateParent, VolumeFromDx, {});
    } else if (updateParent.transformer?.definition.to.includes(PluginStateObject.Format.Cif)) {
      volume = UpdateTarget.apply(updateParent, VolumeFromDensityServerCif, {
        blockHeader: node.params.channel_id || undefined,
      });
    } else {
      console.error(`Unsupported volume format`);
      return undefined;
    }
    return transformAndInstantiateVolume(volume, node);
  },
  volume_representation(
    updateParent: UpdateTarget,
    node: MolstarNode<'volume_representation'>,
    context: MolstarLoadingContext,
  ): UpdateTarget {
    return UpdateTarget.apply(updateParent, VolumeRepresentation3D, {
      ...volumeRepresentationProps(node),
      colorTheme: volumeColorThemeForNode(node, context),
    });
  },
  shape(updateParent: UpdateTarget, node: MolstarSubtree<'shape'>): UpdateTarget | undefined {
    // Same transformers the plugin's own shape format provider uses, so the resulting state
    // object is identical to one produced by opening the file through the UI. MVS only
    // supplies the initial parameter values -- the `volume` node works the same way.
    const transforms = shapeTransforms(node);
    let shape: UpdateTarget;
    let format: ShapeFormat;
    if (updateParent.transformer?.definition.to.includes(PluginStateObject.Format.Vtp)) {
      format = 'vtp';
      shape = UpdateTarget.apply(updateParent, ShapeFromVtp, { transforms });
    } else if (updateParent.transformer?.definition.to.includes(PluginStateObject.Format.Ply)) {
      format = 'ply';
      shape = UpdateTarget.apply(updateParent, ShapeFromPly, { transforms });
    } else if (updateParent.transformer?.definition.to.includes(PluginStateObject.Format.Obj)) {
      format = 'obj';
      shape = UpdateTarget.apply(updateParent, ShapeFromObj, { transforms });
    } else {
      console.error(`Unsupported format for "shape" node`);
      return undefined;
    }
    return UpdateTarget.apply(shape, ShapeRepresentation3D, shapeRepresentationProps(node, format));
  },
  color: undefined, // No action needed, already loaded in `representation`
  color_from_uri: undefined, // No action needed, already loaded in `representation`
  color_from_source: undefined, // No action needed, already loaded in `representation`
  label: undefined, // No action needed, already loaded in `structure`
  label_from_uri(
    updateParent: UpdateTarget,
    node: MolstarNode<'label_from_uri'>,
    context: MolstarLoadingContext,
  ): UpdateTarget {
    const props = labelFromXProps(node, context);
    return UpdateTarget.apply(updateParent, StructureRepresentation3D, props);
  },
  label_from_source(
    updateParent: UpdateTarget,
    node: MolstarNode<'label_from_source'>,
    context: MolstarLoadingContext,
  ): UpdateTarget {
    const props = labelFromXProps(node, context);
    return UpdateTarget.apply(updateParent, StructureRepresentation3D, props);
  },
  focus(updateParent: UpdateTarget, node: MolstarNode<'focus'>, context: MolstarLoadingContext): UpdateTarget {
    context.camera.focuses.push({ target: updateParent.selector, params: node.params });
    return updateParent;
  },
  camera(updateParent: UpdateTarget, node: MolstarNode<'camera'>, context: MolstarLoadingContext): UpdateTarget {
    context.camera.cameraParams = node.params;
    return updateParent;
  },
  canvas(updateParent: UpdateTarget, node: MolstarNode<'canvas'>, context: MolstarLoadingContext): UpdateTarget {
    context.canvas = node;
    return updateParent;
  },
  transition(
    updateParent: UpdateTarget,
    node: MolstarNode<'transition'>,
    context: MolstarLoadingContext,
  ): UpdateTarget {
    context.transition = node;
    return updateParent;
  },
  primitives(
    updateParent: UpdateTarget,
    tree: MolstarSubtree<'primitives'>,
    context: MolstarLoadingContext,
  ): UpdateTarget {
    const refs = getPrimitiveStructureRefs(tree);
    const clip = clippingForNode(tree);
    const data = UpdateTarget.apply(updateParent, MVSInlinePrimitiveData, { node: tree as any });
    UpdateTarget.setMvsDependencies(data, refs); // MVSInlinePrimitiveData must depend on `refs` because it caches positions
    return applyPrimitiveVisuals(data, refs, clip);
  },
  primitives_from_uri(
    updateParent: UpdateTarget,
    tree: MolstarNode<'primitives_from_uri'>,
    context: MolstarLoadingContext,
  ): UpdateTarget {
    const refs = new Set(tree.params.references);
    const clip = clippingForNode(tree);
    const data = UpdateTarget.apply(updateParent, MVSDownloadPrimitiveData, {
      uri: tree.params.uri,
      format: tree.params.format,
    });
    UpdateTarget.setMvsDependencies(data, refs); // MVSInlinePrimitiveData must depend on `refs` because it caches positions
    return applyPrimitiveVisuals(data, refs, clip);
  },
};

function applyPrimitiveVisuals(data: UpdateTarget, refs: Set<string>, clip: Clip.Props | undefined) {
  const mesh = UpdateTarget.setMvsDependencies(
    UpdateTarget.apply(data, MVSBuildPrimitiveShape, { kind: 'mesh', clip }, { state: { isGhost: true } }),
    refs,
  );
  UpdateTarget.apply(mesh, MVSShapeRepresentation3D);
  const labels = UpdateTarget.setMvsDependencies(
    UpdateTarget.apply(data, MVSBuildPrimitiveShape, { kind: 'labels', clip }, { state: { isGhost: true } }),
    refs,
  );
  UpdateTarget.apply(labels, MVSShapeRepresentation3D);
  const lines = UpdateTarget.setMvsDependencies(
    UpdateTarget.apply(data, MVSBuildPrimitiveShape, { kind: 'lines', clip }, { state: { isGhost: true } }),
    refs,
  );
  UpdateTarget.apply(lines, MVSShapeRepresentation3D);
  return data;
}

export type MolstarLoadingExtension<TExtensionContext> = LoadingExtension<
  MolstarTree,
  MolstarLoadingContext,
  TExtensionContext
>;

export const BuiltinLoadingExtensions: MolstarLoadingExtension<any>[] = [
  NonCovalentInteractionsExtension,
  IsHiddenCustomStateExtension,
  VolumeStreamingExtension,
];
