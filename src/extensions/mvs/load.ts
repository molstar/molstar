/**
 * Copyright (c) 2023-2025 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Adam Midlik <midlik@gmail.com>
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Aliaksei Chareshneu <chareshneu.tech@gmail.com>
 */

import { PluginStateSnapshotManager } from '../../mol-plugin-state/manager/snapshots';
import { PluginStateObject } from '../../mol-plugin-state/objects';
import { Download, ParseCif, ParseCcp4 } from '../../mol-plugin-state/transforms/data';
import { CoordinatesFromLammpstraj, CoordinatesFromXtc, CustomModelProperties, CustomStructureProperties, ModelFromTrajectory, StructureComponent, StructureFromModel, TrajectoryFromGRO, TrajectoryFromLammpsTrajData, TrajectoryFromMmCif, TrajectoryFromMOL, TrajectoryFromMOL2, TrajectoryFromPDB, TrajectoryFromSDF, TrajectoryFromXYZ } from '../../mol-plugin-state/transforms/model';
import { StructureRepresentation3D, VolumeRepresentation3D } from '../../mol-plugin-state/transforms/representation';
import { VolumeFromCcp4, VolumeFromDensityServerCif } from '../../mol-plugin-state/transforms/volume';
import { PluginCommands } from '../../mol-plugin/commands';
import { PluginContext } from '../../mol-plugin/context';
import { PluginState } from '../../mol-plugin/state';
import { StateObjectSelector, StateTree } from '../../mol-state';
import { RuntimeContext, Task } from '../../mol-task';
import { MolViewSpec } from './behavior';
import { createPluginStateSnapshotCamera, modifyCanvasProps, resetCanvasProps } from './camera';
import { MVSAnnotationsProvider } from './components/annotation-prop';
import { MVSAnnotationStructureComponent } from './components/annotation-structure-component';
import { MVSAnnotationTooltipsProvider } from './components/annotation-tooltips-prop';
import { CustomLabelProps, CustomLabelRepresentationProvider } from './components/custom-label/representation';
import { CustomTooltipsProvider } from './components/custom-tooltips-prop';
import { IsMVSModelProps, IsMVSModelProvider } from './components/is-mvs-model-prop';
import { getPrimitiveStructureRefs, MVSBuildPrimitiveShape, MVSDownloadPrimitiveData, MVSInlinePrimitiveData, MVSShapeRepresentation3D } from './components/primitives';
import { MVSTrajectoryWithCoordinates } from './components/trajectory';
import { generateStateTransition } from './helpers/animation';
import { IsHiddenCustomStateExtension } from './load-extensions/is-hidden-custom-state';
import { NonCovalentInteractionsExtension } from './load-extensions/non-covalent-interactions';
import { LoadingActions, LoadingExtension, loadTreeVirtual, UpdateTarget } from './load-generic';
import { AnnotationFromSourceKind, AnnotationFromUriKind, collectAnnotationReferences, collectAnnotationTooltips, collectInlineLabels, collectInlineTooltips, colorThemeForNode, componentFromXProps, componentPropsFromSelector, isPhantomComponent, labelFromXProps, makeNearestReprMap, prettyNameFromSelector, representationProps, structureProps, transformAndInstantiateStructure, transformAndInstantiateVolume, volumeColorThemeForNode, volumeRepresentationProps } from './load-helpers';
import { MVSData, MVSData_States, Snapshot, SnapshotMetadata } from './mvs-data';
import { MVSAnimationNode } from './tree/animation/animation-tree';
import { validateTree } from './tree/generic/tree-schema';
import { convertMvsToMolstar, mvsSanityCheck } from './tree/molstar/conversion';
import { MolstarNode, MolstarNodeParams, MolstarSubtree, MolstarTree, MolstarTreeSchema } from './tree/molstar/molstar-tree';
import { MVSTreeSchema } from './tree/mvs/mvs-tree';


export interface MVSLoadOptions {
    /** Add snapshots from MVS into current snapshot list, instead of replacing the list. */
    appendSnapshots?: boolean,
    /** Ignore any camera positioning from the MVS state and keep the current camera position instead, ignore any camera positioning when generating snapshots. */
    keepCamera?: boolean,
    /** Specifies a set of MVS-loading extensions (not a part of standard MVS specification). If undefined, apply all builtin extensions. If `[]`, do not apply builtin extensions. */
    extensions?: MolstarLoadingExtension<any>[],
    /** Run some sanity checks and print potential issues to the console. */
    sanityChecks?: boolean,
    /** Base for resolving relative URLs/URIs. May itself be a relative URL (relative to the window URL). */
    sourceUrl?: string,
    doNotReportErrors?: boolean
}

export function loadMVS(plugin: PluginContext, data: MVSData, options: MVSLoadOptions = {}) {
    const task = Task.create('Load MVS', ctx => _loadMVS(ctx, plugin, data, options));
    return plugin.runTask(task);
}

/** Load a MolViewSpec (MVS) state(s) into the Mol* plugin as plugin state snapshots. */
async function _loadMVS(ctx: RuntimeContext, plugin: PluginContext, data: MVSData, options: MVSLoadOptions = {}) {
    plugin.errorContext.clear('mvs');
    try {
        const mvsExtensionLoaded = plugin.state.hasBehavior(MolViewSpec);
        if (!mvsExtensionLoaded) throw new Error('MolViewSpec extension is not loaded.');

        // Stop any currently running audio
        plugin.managers.markdownExtensions.audio.dispose();

        // Reset canvas props to default so that modifyCanvasProps works as expected
        resetCanvasProps(plugin);

        // console.log(`MVS tree:\n${MVSData.toPrettyString(data)}`)
        const multiData: MVSData_States = data.kind === 'multiple' ? data : MVSData.stateToStates(data);
        const entries: PluginStateSnapshotManager.Entry[] = [];
        for (let i = 0; i < multiData.snapshots.length; i++) {
            const snapshot = multiData.snapshots[i];
            const previousSnapshot = i > 0 ? multiData.snapshots[i - 1] : multiData.snapshots[multiData.snapshots.length - 1];
            validateTree(MVSTreeSchema, snapshot.root, 'MVS');
            if (options.sanityChecks) mvsSanityCheck(snapshot.root);
            const molstarTree = convertMvsToMolstar(snapshot.root, options.sourceUrl);
            validateTree(MolstarTreeSchema, molstarTree, 'Converted Molstar');
            const entry = molstarTreeToEntry(
                plugin,
                molstarTree,
                snapshot.animation,
                { ...snapshot.metadata, previousTransitionDurationMs: previousSnapshot.metadata.transition_duration_ms },
                options
            );
            await assignStateTransition(ctx, plugin, entry, snapshot, options, i, multiData.snapshots.length);
            entries.push(entry);

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
            await PluginCommands.State.Snapshots.Apply(plugin, { id: entries[0].snapshot.id });
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
                    timeoutMs: 10000
                });
            }
        }
        plugin.errorContext.clear('mvs');
    }
}

async function assignStateTransition(ctx: RuntimeContext, plugin: PluginContext, parentEntry: PluginStateSnapshotManager.Entry, parent: Snapshot, options: MVSLoadOptions, snapshotIndex: number, snapshotCount: number) {
    const transitions = await generateStateTransition(ctx, parent, snapshotIndex, snapshotCount);
    if (!transitions?.frames.length) return;

    const animation: PluginState.StateTransition = {
        autoplay: !!transitions.tree.params?.autoplay,
        loop: !!transitions.tree.params?.loop,
        frames: [],
    };

    for (let i = 0; i < transitions.frames.length; i++) {
        const frame = transitions.frames[i];
        const molstarTree = convertMvsToMolstar(frame, options.sourceUrl);
        const entry = molstarTreeToEntry(
            plugin,
            molstarTree,
            parent.animation,
            { ...parent.metadata, previousTransitionDurationMs: transitions.frametimeMs },
            options
        );

        StateTree.reuseTransformParams(entry.snapshot.data!.tree, parentEntry.snapshot.data!.tree);

        animation.frames.push({
            durationInMs: transitions.frametimeMs,
            data: entry.snapshot.data!,
            camera: transitions.tree.params?.include_camera ? entry.snapshot.camera : undefined,
            canvas3d: transitions.tree.params?.include_canvas ? entry.snapshot.canvas3d : undefined,
        });

        if (ctx.shouldUpdate) {
            await ctx.update({ message: `Loading animation for snapshot ${snapshotIndex + 1}/${snapshotCount}...`, current: i + 1, max: transitions.frames.length });
        }
    }

    parentEntry.snapshot.transition = animation;
}

function molstarTreeToEntry(
    plugin: PluginContext,
    tree: MolstarTree,
    animation: MVSAnimationNode<'animation'> | undefined,
    metadata: SnapshotMetadata & { previousTransitionDurationMs?: number },
    options: { keepCamera?: boolean, extensions?: MolstarLoadingExtension<any>[] }
) {
    const context = MolstarLoadingContext.create();
    const snapshot = loadTreeVirtual(plugin, tree, MolstarLoadingActions, context, { replaceExisting: true, extensions: options?.extensions ?? BuiltinLoadingExtensions });
    snapshot.canvas3d = {
        props: plugin.canvas3d ? modifyCanvasProps(plugin.canvas3d.props, context.canvas, animation) : undefined,
    };
    if (!options?.keepCamera) {
        snapshot.camera = createPluginStateSnapshotCamera(plugin, context, metadata);
    }
    snapshot.durationInMs = metadata.linger_duration_ms + (metadata.previousTransitionDurationMs ?? 0);

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
    annotationMap: Map<MolstarNode<AnnotationFromUriKind | AnnotationFromSourceKind>, string>,
    /** Maps each node (on 'structure' or lower level) to its nearest 'representation' node */
    nearestReprMap?: Map<MolstarNode, MolstarNode<'representation'>>,
    camera: {
        cameraParams?: MolstarNodeParams<'camera'>,
        focuses: { target: StateObjectSelector, params: MolstarNodeParams<'focus'> }[],
    },
    canvas?: MolstarNode<'canvas'>,
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
                return updateParent;
            case 'map':
                return UpdateTarget.apply(updateParent, ParseCcp4, {});
            default:
                console.error(`Unknown format in "parse" node: "${format}"`);
                return undefined;
        }
    },
    coordinates(updateParent: UpdateTarget, node: MolstarNode<'coordinates'>): UpdateTarget | undefined {
        const format = node.params.format;
        switch (format) {
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
                return UpdateTarget.apply(updateParent, TrajectoryFromPDB, { isPdbqt: format === 'pdbqt' });
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
    trajectory_with_coordinates(updateParent: UpdateTarget, node: MolstarNode<'trajectory_with_coordinates'>): UpdateTarget | undefined {
        const result = UpdateTarget.apply(updateParent, MVSTrajectoryWithCoordinates, {
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
            autoAttach: [
                IsMVSModelProvider.descriptor.name,
                MVSAnnotationsProvider.descriptor.name,
            ],
        });
        return model;
    },
    structure(updateParent: UpdateTarget, node: MolstarSubtree<'structure'>, context: MolstarLoadingContext): UpdateTarget {
        const props = structureProps(node);
        const struct = UpdateTarget.apply(updateParent, StructureFromModel, props);
        const transformed = transformAndInstantiateStructure(struct, node);
        const annotationTooltips = collectAnnotationTooltips(node, context);
        const inlineTooltips = collectInlineTooltips(node, context);
        if (annotationTooltips.length + inlineTooltips.length > 0) {
            UpdateTarget.apply(struct, CustomStructureProperties, {
                properties: {
                    [MVSAnnotationTooltipsProvider.descriptor.name]: { tooltips: annotationTooltips },
                    [CustomTooltipsProvider.descriptor.name]: { tooltips: inlineTooltips },
                },
                autoAttach: [
                    MVSAnnotationTooltipsProvider.descriptor.name,
                    CustomTooltipsProvider.descriptor.name,
                ],
            });
        }
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
        return transformAndInstantiateStructure(UpdateTarget.apply(updateParent, StructureComponent, {
            type: componentPropsFromSelector(selector),
            label: prettyNameFromSelector(selector),
            nullIfEmpty: false,
        }), node);
    },
    component_from_uri(updateParent: UpdateTarget, node: MolstarSubtree<'component_from_uri'>, context: MolstarLoadingContext): UpdateTarget | undefined {
        if (isPhantomComponent(node)) return undefined;
        const props = componentFromXProps(node, context);
        return transformAndInstantiateStructure(UpdateTarget.apply(updateParent, MVSAnnotationStructureComponent, props), node);
    },
    component_from_source(updateParent: UpdateTarget, node: MolstarSubtree<'component_from_source'>, context: MolstarLoadingContext): UpdateTarget | undefined {
        if (isPhantomComponent(node)) return undefined;
        const props = componentFromXProps(node, context);
        return transformAndInstantiateStructure(UpdateTarget.apply(updateParent, MVSAnnotationStructureComponent, props), node);
    },
    representation(updateParent: UpdateTarget, node: MolstarNode<'representation'>, context: MolstarLoadingContext): UpdateTarget {
        return UpdateTarget.apply(updateParent, StructureRepresentation3D, {
            ...representationProps(node),
            colorTheme: colorThemeForNode(node, context),
        });
    },
    volume(updateParent: UpdateTarget, node: MolstarNode<'volume'>): UpdateTarget | undefined {
        let volume: UpdateTarget;
        if (updateParent.transformer?.definition.to.includes(PluginStateObject.Format.Ccp4)) {
            volume = UpdateTarget.apply(updateParent, VolumeFromCcp4, {});
        } else if (updateParent.transformer?.definition.to.includes(PluginStateObject.Format.Cif)) {
            volume = UpdateTarget.apply(updateParent, VolumeFromDensityServerCif, { blockHeader: node.params.channel_id || undefined });
        } else {
            console.error(`Unsupported volume format`);
            return undefined;
        }
        return transformAndInstantiateVolume(volume, node);
    },
    volume_representation(updateParent: UpdateTarget, node: MolstarNode<'volume_representation'>, context: MolstarLoadingContext): UpdateTarget {
        return UpdateTarget.apply(updateParent, VolumeRepresentation3D, {
            ...volumeRepresentationProps(node),
            colorTheme: volumeColorThemeForNode(node, context),
        });
    },
    color: undefined, // No action needed, already loaded in `representation`
    color_from_uri: undefined, // No action needed, already loaded in `representation`
    color_from_source: undefined, // No action needed, already loaded in `representation`
    label: undefined, // No action needed, already loaded in `structure`
    label_from_uri(updateParent: UpdateTarget, node: MolstarNode<'label_from_uri'>, context: MolstarLoadingContext): UpdateTarget {
        const props = labelFromXProps(node, context);
        return UpdateTarget.apply(updateParent, StructureRepresentation3D, props);
    },
    label_from_source(updateParent: UpdateTarget, node: MolstarNode<'label_from_source'>, context: MolstarLoadingContext): UpdateTarget {
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
    primitives(updateParent: UpdateTarget, tree: MolstarSubtree<'primitives'>, context: MolstarLoadingContext): UpdateTarget {
        const refs = getPrimitiveStructureRefs(tree);
        const data = UpdateTarget.apply(updateParent, MVSInlinePrimitiveData, { node: tree as any });
        return applyPrimitiveVisuals(data, refs);
    },
    primitives_from_uri(updateParent: UpdateTarget, tree: MolstarNode<'primitives_from_uri'>, context: MolstarLoadingContext): UpdateTarget {
        const data = UpdateTarget.apply(updateParent, MVSDownloadPrimitiveData, { uri: tree.params.uri, format: tree.params.format });
        return applyPrimitiveVisuals(data, new Set(tree.params.references));
    },
};

function applyPrimitiveVisuals(data: UpdateTarget, refs: Set<string>) {
    const mesh = UpdateTarget.setMvsDependencies(UpdateTarget.apply(data, MVSBuildPrimitiveShape, { kind: 'mesh' }, { state: { isGhost: true } }), refs);
    UpdateTarget.apply(mesh, MVSShapeRepresentation3D);
    const labels = UpdateTarget.setMvsDependencies(UpdateTarget.apply(data, MVSBuildPrimitiveShape, { kind: 'labels' }, { state: { isGhost: true } }), refs);
    UpdateTarget.apply(labels, MVSShapeRepresentation3D);
    const lines = UpdateTarget.setMvsDependencies(UpdateTarget.apply(data, MVSBuildPrimitiveShape, { kind: 'lines' }, { state: { isGhost: true } }), refs);
    UpdateTarget.apply(lines, MVSShapeRepresentation3D);
    return data;
}

export type MolstarLoadingExtension<TExtensionContext> = LoadingExtension<MolstarTree, MolstarLoadingContext, TExtensionContext>;

export const BuiltinLoadingExtensions: MolstarLoadingExtension<any>[] = [
    NonCovalentInteractionsExtension,
    IsHiddenCustomStateExtension,
];
