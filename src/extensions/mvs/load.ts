/**
 * Copyright (c) 2023 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Adam Midlik <midlik@gmail.com>
 */

import { Download, ParseCif } from '../../mol-plugin-state/transforms/data';
import { CustomModelProperties, CustomStructureProperties, ModelFromTrajectory, StructureComponent, StructureFromModel, TrajectoryFromMmCif, TrajectoryFromPDB, TransformStructureConformation } from '../../mol-plugin-state/transforms/model';
import { StructureRepresentation3D } from '../../mol-plugin-state/transforms/representation';
import { PluginContext } from '../../mol-plugin/context';
import { StateBuilder, StateObjectSelector } from '../../mol-state';
import { AnnotationsProvider } from './additions/annotation-prop';
import { AnnotationStructureComponent } from './additions/annotation-structure-component';
import { AnnotationTooltipsProvider } from './additions/annotation-tooltips-prop';
import { CustomLabelProps, CustomLabelRepresentationProvider } from './additions/custom-label/representation';
import { CustomTooltipsProvider } from './additions/custom-tooltips-prop';
import { setCamera, setCanvas, setFocus } from './camera';
import { canonicalJsonString } from './helpers/utils';
import { AnnotationFromSourceKind, AnnotationFromUriKind, LoadingActions, collectAnnotationReferences, collectAnnotationTooltips, collectInlineTooltips, colorThemeForNode, componentFromXProps, componentPropsFromSelector, isPhantomComponent, labelFromXProps, loadTree, makeNearestReprMap, representationProps, structureProps, transformProps } from './load-helpers';
import { MVSData } from './mvs-data';
import { ParamsOfKind, SubTreeOfKind, validateTree } from './tree/generic/tree-schema';
import { convertMvsToMolstar } from './tree/molstar/conversion';
import { MolstarNode, MolstarTree, MolstarTreeSchema } from './tree/molstar/molstar-tree';
import { MVSTreeSchema } from './tree/mvs/mvs-tree';


/** Load a MolViewSpec (MVS) tree into the Mol* plugin.
 * If `deletePrevious`, remove all objects in the current Mol* state; otherwise add to the current state.
 */
export async function loadMVS(plugin: PluginContext, data: MVSData | string, deletePrevious: boolean) {
    if (typeof data === 'string') {
        data = JSON.parse(data) as MVSData;
    }
    // console.log(`MVS tree (v${data.version}):\n${treeToString(data.root)}`);
    validateTree(MVSTreeSchema, data.root, 'MVS');
    const molstarTree = convertMvsToMolstar(data.root);
    // console.log(`Converted MolStar tree:\n${treeToString(molstarTree)}`);
    validateTree(MolstarTreeSchema, molstarTree, 'Converted Molstar');
    await loadMolstarTree(plugin, molstarTree, deletePrevious);
}


/** Load a `MolstarTree` into the Mol* plugin.
 * If `deletePrevious`, remove all objects in the current Mol* state; otherwise add to the current state.
 */
async function loadMolstarTree(plugin: PluginContext, tree: MolstarTree, deletePrevious: boolean) {
    const context: MolstarLoadingContext = {};

    await loadTree(plugin, tree, MolstarLoadingActions, context, deletePrevious);

    setCanvas(plugin, context.canvas);
    if (context.focus?.kind === 'camera') {
        await setCamera(plugin, context.focus.params);
    } else if (context.focus?.kind === 'focus') {
        await setFocus(plugin, context.focus.focusTarget, context.focus.params);
    } else {
        await setFocus(plugin, undefined, undefined);
    }
}

/** Mutable context for loading a `MolstarTree`, available throughout the loading. */
export interface MolstarLoadingContext {
    /** Maps `*_from_[uri|source]` nodes to annotationId they should reference */
    annotationMap?: Map<MolstarNode<AnnotationFromUriKind | AnnotationFromSourceKind>, string>,
    /** Maps each node (on 'structure' or lower level) to its nearest 'representation' node */
    nearestReprMap?: Map<MolstarNode, MolstarNode<'representation'>>,
    focus?: { kind: 'camera', params: ParamsOfKind<MolstarTree, 'camera'> } | { kind: 'focus', focusTarget: StateObjectSelector, params: ParamsOfKind<MolstarTree, 'focus'> },
    canvas?: ParamsOfKind<MolstarTree, 'canvas'>,
}


/** Loading actions for loading a `MolstarTree`, per node kind. */
const MolstarLoadingActions: LoadingActions<MolstarTree, MolstarLoadingContext> = {
    root(update: StateBuilder.Root, msParent: StateObjectSelector, node: MolstarNode<'root'>, context: MolstarLoadingContext): StateObjectSelector {
        context.nearestReprMap = makeNearestReprMap(node);
        return msParent;
    },
    download(update: StateBuilder.Root, msParent: StateObjectSelector, node: MolstarNode<'download'>): StateObjectSelector {
        return update.to(msParent).apply(Download, {
            url: node.params.url,
            isBinary: node.params.is_binary,
        }).selector;
    },
    parse(update: StateBuilder.Root, msParent: StateObjectSelector, node: MolstarNode<'parse'>): StateObjectSelector | undefined {
        const format = node.params.format;
        if (format === 'cif') {
            return update.to(msParent).apply(ParseCif, {}).selector;
        } else if (format === 'pdb') {
            return msParent;
        } else {
            console.error(`Unknown format in "parse" node: "${format}"`);
            return undefined;
        }
    },
    trajectory(update: StateBuilder.Root, msParent: StateObjectSelector, node: MolstarNode<'trajectory'>): StateObjectSelector | undefined {
        const format = node.params.format;
        if (format === 'cif') {
            return update.to(msParent).apply(TrajectoryFromMmCif, {
                blockHeader: node.params.block_header ?? '', // Must set to '' because just undefined would get overwritten by createDefaults
                blockIndex: node.params.block_index ?? undefined,
            }).selector;
        } else if (format === 'pdb') {
            return update.to(msParent).apply(TrajectoryFromPDB, {}).selector;
        } else {
            console.error(`Unknown format in "trajectory" node: "${format}"`);
            return undefined;
        }
    },
    model(update: StateBuilder.Root, msParent: StateObjectSelector, node: SubTreeOfKind<MolstarTree, 'model'>, context: MolstarLoadingContext): StateObjectSelector {
        const annotations = collectAnnotationReferences(node, context);
        return update.to(msParent)
            .apply(ModelFromTrajectory, {
                modelIndex: node.params.model_index,
            })
            .apply(CustomModelProperties, {
                properties: {
                    [AnnotationsProvider.descriptor.name]: { annotations }
                },
                autoAttach: [
                    AnnotationsProvider.descriptor.name
                ],
            }).selector;
    },
    structure(update: StateBuilder.Root, msParent: StateObjectSelector, node: SubTreeOfKind<MolstarTree, 'structure'>, context: MolstarLoadingContext): StateObjectSelector {
        const props = structureProps(node);
        let result: StateObjectSelector = update.to(msParent).apply(StructureFromModel, props).selector;
        for (const t of transformProps(node)) {
            result = update.to(result).apply(TransformStructureConformation, t).selector;
        }
        const annotationTooltips = collectAnnotationTooltips(node, context);
        const inlineTooltips = collectInlineTooltips(node, context);
        if (annotationTooltips.length + inlineTooltips.length > 0) {
            update.to(result).apply(CustomStructureProperties, {
                properties: {
                    [AnnotationTooltipsProvider.descriptor.name]: { tooltips: annotationTooltips },
                    [CustomTooltipsProvider.descriptor.name]: { tooltips: inlineTooltips },
                },
                autoAttach: [
                    AnnotationTooltipsProvider.descriptor.name,
                    CustomTooltipsProvider.descriptor.name,
                ],
            });
        }
        return result;
    },
    tooltip: undefined, // No action needed, already loaded in `structure`
    tooltip_from_uri: undefined, // No action needed, already loaded in `structure`
    tooltip_from_source: undefined, // No action needed, already loaded in `structure`
    component(update: StateBuilder.Root, msParent: StateObjectSelector, node: SubTreeOfKind<MolstarTree, 'component'>): StateObjectSelector | undefined {
        if (isPhantomComponent(node)) {
            return msParent;
        }
        const selector = node.params.selector;
        return update.to(msParent).apply(StructureComponent, {
            type: componentPropsFromSelector(selector),
            label: canonicalJsonString(selector),
            nullIfEmpty: false,
        }).selector;
    },
    component_from_uri(update: StateBuilder.Root, msParent: StateObjectSelector, node: SubTreeOfKind<MolstarTree, 'component_from_uri'>, context: MolstarLoadingContext): StateObjectSelector | undefined {
        if (isPhantomComponent(node)) return undefined;
        const props = componentFromXProps(node, context);
        return update.to(msParent).apply(AnnotationStructureComponent, props).selector;
    },
    component_from_source(update: StateBuilder.Root, msParent: StateObjectSelector, node: SubTreeOfKind<MolstarTree, 'component_from_source'>, context: MolstarLoadingContext): StateObjectSelector | undefined {
        if (isPhantomComponent(node)) return undefined;
        const props = componentFromXProps(node, context);
        return update.to(msParent).apply(AnnotationStructureComponent, props).selector;
    },
    representation(update: StateBuilder.Root, msParent: StateObjectSelector, node: MolstarNode<'representation'>, context: MolstarLoadingContext): StateObjectSelector {
        return update.to(msParent).apply(StructureRepresentation3D, {
            ...representationProps(node.params),
            colorTheme: colorThemeForNode(node, context),
        }).selector;
    },
    color: undefined, // No action needed, already loaded in `structure`
    color_from_uri: undefined, // No action needed, already loaded in `structure`
    color_from_source: undefined, // No action needed, already loaded in `structure`
    label(update: StateBuilder.Root, msParent: StateObjectSelector, node: MolstarNode<'label'>, context: MolstarLoadingContext): StateObjectSelector {
        const item: CustomLabelProps['items'][number] = {
            text: node.params.text,
            position: { name: 'selection', params: {} },
        };
        const nearestReprNode = context.nearestReprMap?.get(node);
        return update.to(msParent).apply(StructureRepresentation3D, {
            type: {
                name: CustomLabelRepresentationProvider.name,
                params: { items: [item] } satisfies Partial<CustomLabelProps>
            },
            colorTheme: colorThemeForNode(nearestReprNode, context),
        }).selector;
    },
    label_from_uri(update: StateBuilder.Root, msParent: StateObjectSelector, node: MolstarNode<'label_from_uri'>, context: MolstarLoadingContext): StateObjectSelector {
        const props = labelFromXProps(node, context);
        return update.to(msParent).apply(StructureRepresentation3D, props).selector;
    },
    label_from_source(update: StateBuilder.Root, msParent: StateObjectSelector, node: MolstarNode<'label_from_source'>, context: MolstarLoadingContext): StateObjectSelector {
        const props = labelFromXProps(node, context);
        return update.to(msParent).apply(StructureRepresentation3D, props).selector;
    },
    focus(update: StateBuilder.Root, msParent: StateObjectSelector, node: MolstarNode<'focus'>, context: MolstarLoadingContext): StateObjectSelector {
        context.focus = { kind: 'focus', focusTarget: msParent, params: node.params };
        return msParent;
    },
    camera(update: StateBuilder.Root, msParent: StateObjectSelector, node: MolstarNode<'camera'>, context: MolstarLoadingContext): StateObjectSelector {
        context.focus = { kind: 'camera', params: node.params };
        return msParent;
    },
    canvas(update: StateBuilder.Root, msParent: StateObjectSelector, node: MolstarNode<'canvas'>, context: MolstarLoadingContext): StateObjectSelector {
        context.canvas = node.params;
        return msParent;
    },
};
