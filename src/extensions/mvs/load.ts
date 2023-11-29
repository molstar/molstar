/**
 * Copyright (c) 2023 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Adam Midlik <midlik@gmail.com>
 */

import { Download, ParseCif } from '../../mol-plugin-state/transforms/data';
import { CustomModelProperties, CustomStructureProperties, ModelFromTrajectory, StructureComponent, StructureFromModel, TrajectoryFromMmCif, TrajectoryFromPDB, TransformStructureConformation } from '../../mol-plugin-state/transforms/model';
import { StructureRepresentation3D } from '../../mol-plugin-state/transforms/representation';
import { PluginContext } from '../../mol-plugin/context';
import { StateObjectSelector } from '../../mol-state';
import { canonicalJsonString } from '../../mol-util/json';
import { MolViewSpec } from './behavior';
import { setCamera, setCanvas, setFocus } from './camera';
import { MVSAnnotationsProvider } from './components/annotation-prop';
import { MVSAnnotationStructureComponent } from './components/annotation-structure-component';
import { MVSAnnotationTooltipsProvider } from './components/annotation-tooltips-prop';
import { CustomLabelProps, CustomLabelRepresentationProvider } from './components/custom-label/representation';
import { CustomTooltipsProvider } from './components/custom-tooltips-prop';
import { AnnotationFromSourceKind, AnnotationFromUriKind, LoadingActions, UpdateTarget, collectAnnotationReferences, collectAnnotationTooltips, collectInlineTooltips, colorThemeForNode, componentFromXProps, componentPropsFromSelector, isPhantomComponent, labelFromXProps, loadTree, makeNearestReprMap, representationProps, structureProps, transformProps } from './load-helpers';
import { MVSData } from './mvs-data';
import { ParamsOfKind, SubTreeOfKind, validateTree } from './tree/generic/tree-schema';
import { convertMvsToMolstar, mvsSanityCheck } from './tree/molstar/conversion';
import { MolstarNode, MolstarTree, MolstarTreeSchema } from './tree/molstar/molstar-tree';
import { MVSTreeSchema } from './tree/mvs/mvs-tree';


/** Load a MolViewSpec (MVS) tree into the Mol* plugin.
 * If `options.replaceExisting`, remove all objects in the current Mol* state; otherwise add to the current state.
 * If `options.sanityChecks`, run some sanity checks and print potential issues to the console. */
export async function loadMVS(plugin: PluginContext, data: MVSData, options: { replaceExisting?: boolean, sanityChecks?: boolean } = {}) {
    try {
        // console.log(`MVS tree (v${data.version}):\n${treeToString(data.root)}`);
        validateTree(MVSTreeSchema, data.root, 'MVS');
        if (options.sanityChecks) mvsSanityCheck(data.root);
        const molstarTree = convertMvsToMolstar(data.root);
        // console.log(`Converted MolStar tree:\n${treeToString(molstarTree)}`);
        validateTree(MolstarTreeSchema, molstarTree, 'Converted Molstar');
        await loadMolstarTree(plugin, molstarTree, options);
    } catch (err) {
        plugin.log.error(`${err}`);
        throw err;
    }
}


/** Load a `MolstarTree` into the Mol* plugin.
 * If `replaceExisting`, remove all objects in the current Mol* state; otherwise add to the current state. */
async function loadMolstarTree(plugin: PluginContext, tree: MolstarTree, options?: { replaceExisting?: boolean }) {
    const mvsExtensionLoaded = plugin.state.hasBehavior(MolViewSpec);
    if (!mvsExtensionLoaded) throw new Error('MolViewSpec extension is not loaded.');

    const context: MolstarLoadingContext = {};

    await loadTree(plugin, tree, MolstarLoadingActions, context, options);

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
        if (format === 'cif') {
            return UpdateTarget.apply(updateParent, ParseCif, {});
        } else if (format === 'pdb') {
            return updateParent;
        } else {
            console.error(`Unknown format in "parse" node: "${format}"`);
            return undefined;
        }
    },
    trajectory(updateParent: UpdateTarget, node: MolstarNode<'trajectory'>): UpdateTarget | undefined {
        const format = node.params.format;
        if (format === 'cif') {
            return UpdateTarget.apply(updateParent, TrajectoryFromMmCif, {
                blockHeader: node.params.block_header ?? '', // Must set to '' because just undefined would get overwritten by createDefaults
                blockIndex: node.params.block_index ?? undefined,
            });
        } else if (format === 'pdb') {
            return UpdateTarget.apply(updateParent, TrajectoryFromPDB, {});
        } else {
            console.error(`Unknown format in "trajectory" node: "${format}"`);
            return undefined;
        }
    },
    model(updateParent: UpdateTarget, node: SubTreeOfKind<MolstarTree, 'model'>, context: MolstarLoadingContext): UpdateTarget {
        const annotations = collectAnnotationReferences(node, context);
        const model = UpdateTarget.apply(updateParent, ModelFromTrajectory, {
            modelIndex: node.params.model_index,
        });
        UpdateTarget.apply(model, CustomModelProperties, {
            properties: {
                [MVSAnnotationsProvider.descriptor.name]: { annotations }
            },
            autoAttach: [
                MVSAnnotationsProvider.descriptor.name
            ],
        });
        return model;
    },
    structure(updateParent: UpdateTarget, node: SubTreeOfKind<MolstarTree, 'structure'>, context: MolstarLoadingContext): UpdateTarget {
        const props = structureProps(node);
        const struct = UpdateTarget.apply(updateParent, StructureFromModel, props);
        let transformed = struct;
        for (const t of transformProps(node)) {
            transformed = UpdateTarget.apply(transformed, TransformStructureConformation, t); // applying to the result of previous transform, to get the correct transform order
        }
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
        return struct;
    },
    tooltip: undefined, // No action needed, already loaded in `structure`
    tooltip_from_uri: undefined, // No action needed, already loaded in `structure`
    tooltip_from_source: undefined, // No action needed, already loaded in `structure`
    component(updateParent: UpdateTarget, node: SubTreeOfKind<MolstarTree, 'component'>): UpdateTarget | undefined {
        if (isPhantomComponent(node)) {
            return updateParent;
        }
        const selector = node.params.selector;
        return UpdateTarget.apply(updateParent, StructureComponent, {
            type: componentPropsFromSelector(selector),
            label: canonicalJsonString(selector),
            nullIfEmpty: false,
        });
    },
    component_from_uri(updateParent: UpdateTarget, node: SubTreeOfKind<MolstarTree, 'component_from_uri'>, context: MolstarLoadingContext): UpdateTarget | undefined {
        if (isPhantomComponent(node)) return undefined;
        const props = componentFromXProps(node, context);
        return UpdateTarget.apply(updateParent, MVSAnnotationStructureComponent, props);
    },
    component_from_source(updateParent: UpdateTarget, node: SubTreeOfKind<MolstarTree, 'component_from_source'>, context: MolstarLoadingContext): UpdateTarget | undefined {
        if (isPhantomComponent(node)) return undefined;
        const props = componentFromXProps(node, context);
        return UpdateTarget.apply(updateParent, MVSAnnotationStructureComponent, props);
    },
    representation(updateParent: UpdateTarget, node: MolstarNode<'representation'>, context: MolstarLoadingContext): UpdateTarget {
        return UpdateTarget.apply(updateParent, StructureRepresentation3D, {
            ...representationProps(node.params),
            colorTheme: colorThemeForNode(node, context),
        });
    },
    color: undefined, // No action needed, already loaded in `structure`
    color_from_uri: undefined, // No action needed, already loaded in `structure`
    color_from_source: undefined, // No action needed, already loaded in `structure`
    label(updateParent: UpdateTarget, node: MolstarNode<'label'>, context: MolstarLoadingContext): UpdateTarget {
        const item: CustomLabelProps['items'][number] = {
            text: node.params.text,
            position: { name: 'selection', params: {} },
        };
        const nearestReprNode = context.nearestReprMap?.get(node);
        return UpdateTarget.apply(updateParent, StructureRepresentation3D, {
            type: {
                name: CustomLabelRepresentationProvider.name,
                params: { items: [item] } satisfies Partial<CustomLabelProps>
            },
            colorTheme: colorThemeForNode(nearestReprNode, context),
        });
    },
    label_from_uri(updateParent: UpdateTarget, node: MolstarNode<'label_from_uri'>, context: MolstarLoadingContext): UpdateTarget {
        const props = labelFromXProps(node, context);
        return UpdateTarget.apply(updateParent, StructureRepresentation3D, props);
    },
    label_from_source(updateParent: UpdateTarget, node: MolstarNode<'label_from_source'>, context: MolstarLoadingContext): UpdateTarget {
        const props = labelFromXProps(node, context);
        return UpdateTarget.apply(updateParent, StructureRepresentation3D, props);
    },
    focus(updateParent: UpdateTarget, node: MolstarNode<'focus'>, context: MolstarLoadingContext): UpdateTarget {
        context.focus = { kind: 'focus', focusTarget: updateParent.selector, params: node.params };
        return updateParent;
    },
    camera(updateParent: UpdateTarget, node: MolstarNode<'camera'>, context: MolstarLoadingContext): UpdateTarget {
        context.focus = { kind: 'camera', params: node.params };
        return updateParent;
    },
    canvas(updateParent: UpdateTarget, node: MolstarNode<'canvas'>, context: MolstarLoadingContext): UpdateTarget {
        context.canvas = node.params;
        return updateParent;
    },
};
