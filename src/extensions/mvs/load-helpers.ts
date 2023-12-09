/**
 * Copyright (c) 2023 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Adam Midlik <midlik@gmail.com>
 */

import { Mat3, Mat4, Vec3 } from '../../mol-math/linear-algebra';
import { StructureComponentParams } from '../../mol-plugin-state/helpers/structure-component';
import { StructureFromModel, TransformStructureConformation } from '../../mol-plugin-state/transforms/model';
import { StructureRepresentation3D } from '../../mol-plugin-state/transforms/representation';
import { PluginContext } from '../../mol-plugin/context';
import { StateBuilder, StateObject, StateObjectSelector, StateTransform, StateTransformer } from '../../mol-state';
import { arrayDistinct } from '../../mol-util/array';
import { canonicalJsonString } from '../../mol-util/json';
import { MVSAnnotationColorThemeProps, MVSAnnotationColorThemeProvider } from './components/annotation-color-theme';
import { MVSAnnotationLabelRepresentationProvider } from './components/annotation-label/representation';
import { MVSAnnotationSpec } from './components/annotation-prop';
import { MVSAnnotationStructureComponentProps } from './components/annotation-structure-component';
import { MVSAnnotationTooltipsProps } from './components/annotation-tooltips-prop';
import { CustomTooltipsProps } from './components/custom-tooltips-prop';
import { MultilayerColorThemeName, MultilayerColorThemeProps, NoColor } from './components/multilayer-color-theme';
import { SelectorAll } from './components/selector';
import { rowToExpression, rowsToExpression } from './helpers/selections';
import { ElementOfSet, decodeColor, isDefined, stringHash } from './helpers/utils';
import { MolstarLoadingContext } from './load';
import { Kind, ParamsOfKind, SubTree, SubTreeOfKind, Tree, getChildren } from './tree/generic/tree-schema';
import { dfs } from './tree/generic/tree-utils';
import { MolstarKind, MolstarNode, MolstarTree } from './tree/molstar/molstar-tree';
import { DefaultColor } from './tree/mvs/mvs-defaults';


/** Function responsible for loading a tree node `node` into Mol*.
 * Should apply changes within `updateParent.update` but not commit them.
 * Should modify `context` accordingly, if it is needed for loading other nodes later.
 * `updateParent.selector` is the result of loading the node's parent into Mol* state hierarchy (or the hierarchy root in case of root node). */
export type LoadingAction<TNode extends Tree, TContext> = (updateParent: UpdateTarget, node: TNode, context: TContext) => UpdateTarget | undefined

/** Loading actions for loading a tree into Mol*, per node kind. */
export type LoadingActions<TTree extends Tree, TContext> = { [kind in Kind<SubTree<TTree>>]?: LoadingAction<SubTreeOfKind<TTree, kind>, TContext> }

/** Load a tree into Mol*, by applying loading actions in DFS order and then commiting at once.
 * If `options.replaceExisting`, remove all objects in the current Mol* state; otherwise add to the current state. */
export async function loadTree<TTree extends Tree, TContext>(plugin: PluginContext, tree: TTree, loadingActions: LoadingActions<TTree, TContext>, context: TContext, options?: { replaceExisting?: boolean }) {
    const mapping = new Map<SubTree<TTree>, UpdateTarget | undefined>();
    const updateRoot: UpdateTarget = UpdateTarget.create(plugin, options?.replaceExisting ?? false);
    if (options?.replaceExisting) {
        UpdateTarget.deleteChildren(updateRoot);
    }
    dfs<TTree>(tree, (node, parent) => {
        const kind: Kind<typeof node> = node.kind;
        const action = loadingActions[kind] as LoadingAction<typeof node, TContext> | undefined;
        if (action) {
            const updateParent = parent ? mapping.get(parent) : updateRoot;
            if (updateParent) {
                const msNode = action(updateParent, node, context);
                mapping.set(node, msNode);
            } else {
                console.warn(`No target found for this "${node.kind}" node`);
                return;
            }
        }
    });
    await UpdateTarget.commit(updateRoot);
}


/** A wrapper for updating Mol* state, while using deterministic transform refs.
 * ```
 * updateTarget = UpdateTarget.create(plugin); // like update = plugin.build();
 * UpdateTarget.apply(updateTarget, transformer, params); // like update.to(selector).apply(transformer, params);
 * await UpdateTarget.commit(updateTarget); // like await update.commit();
 * ```
 */
export interface UpdateTarget {
    readonly update: StateBuilder.Root,
    readonly selector: StateObjectSelector,
    readonly refManager: RefManager,
}
export const UpdateTarget = {
    /** Create a new update, with `selector` pointing to the root. */
    create(plugin: PluginContext, replaceExisting: boolean): UpdateTarget {
        const update = plugin.build();
        const msTarget = update.toRoot().selector;
        const refManager = new RefManager(plugin, replaceExisting);
        return { update, selector: msTarget, refManager };
    },
    /** Add a child node to `target.selector`, return a new `UpdateTarget` pointing to the new child. */
    apply<A extends StateObject, B extends StateObject, P extends {}>(target: UpdateTarget, transformer: StateTransformer<A, B, P>, params?: Partial<P>, options?: Partial<StateTransform.Options>): UpdateTarget {
        let refSuffix: string = transformer.id;
        if (transformer.id === StructureRepresentation3D.id) {
            const reprType = (params as any)?.type?.name ?? '';
            refSuffix += `:${reprType}`;
        }
        const ref = target.refManager.getChildRef(target.selector, refSuffix);
        const msResult = target.update.to(target.selector).apply(transformer, params, { ...options, ref }).selector;
        return { ...target, selector: msResult };
    },
    /** Delete all children of `target.selector`. */
    deleteChildren(target: UpdateTarget): UpdateTarget {
        const children = target.update.currentTree.children.get(target.selector.ref);
        children.forEach(child => target.update.delete(child));
        return target;
    },
    /** Commit all changes done in the current update. */
    commit(target: UpdateTarget): Promise<void> {
        return target.update.commit();
    },
};

/** Manages transform refs in a deterministic way. Uses refs like !mvs:3ce3664304d32c5d:0 */
class RefManager {
    /** For each hash (e.g. 3ce3664304d32c5d), store the number of already used refs with that hash. */
    private _counter: Record<string, number> = {};
    constructor(plugin: PluginContext, replaceExisting: boolean) {
        if (!replaceExisting) {
            plugin.state.data.cells.forEach(cell => {
                const ref = cell.transform.ref;
                if (ref.startsWith('!mvs:')) {
                    const [_, hash, idNumber] = ref.split(':');
                    const nextIdNumber = parseInt(idNumber) + 1;
                    if (nextIdNumber > (this._counter[hash] ?? 0)) {
                        this._counter[hash] = nextIdNumber;
                    }
                }
            });
        }
    }
    /** Return ref for a new node with given `hash`; update the counter accordingly. */
    private nextRef(hash: string): string {
        this._counter[hash] ??= 0;
        const idNumber = this._counter[hash]++;
        return `!mvs:${hash}:${idNumber}`;
    }
    /** Return ref for a new node based on parent and desired suffix. */
    getChildRef(parent: StateObjectSelector, suffix: string): string {
        const hashBase = parent.ref.replace(/^!mvs:/, '') + ':' + suffix;
        const hash = stringHash(hashBase);
        const result = this.nextRef(hash);
        return result;
    }
}


export const AnnotationFromUriKinds = new Set(['color_from_uri', 'component_from_uri', 'label_from_uri', 'tooltip_from_uri'] satisfies MolstarKind[]);
export type AnnotationFromUriKind = ElementOfSet<typeof AnnotationFromUriKinds>

export const AnnotationFromSourceKinds = new Set(['color_from_source', 'component_from_source', 'label_from_source', 'tooltip_from_source'] satisfies MolstarKind[]);
export type AnnotationFromSourceKind = ElementOfSet<typeof AnnotationFromSourceKinds>


/** Return a 4x4 matrix representing a rotation followed by a translation */
export function transformFromRotationTranslation(rotation: number[] | null | undefined, translation: number[] | null | undefined): Mat4 {
    if (rotation && rotation.length !== 9) throw new Error(`'rotation' param for 'transform' node must be array of 9 elements, found ${rotation}`);
    if (translation && translation.length !== 3) throw new Error(`'translation' param for 'transform' node must be array of 3 elements, found ${translation}`);
    const T = Mat4.identity();
    if (rotation) {
        const rotMatrix = Mat3.fromArray(Mat3(), rotation, 0);
        ensureRotationMatrix(rotMatrix, rotMatrix);
        Mat4.fromMat3(T, rotMatrix);
    }
    if (translation) {
        Mat4.setTranslation(T, Vec3.fromArray(Vec3(), translation, 0));
    }
    if (!Mat4.isRotationAndTranslation(T)) throw new Error(`'rotation' param for 'transform' is not a valid rotation matrix: ${rotation}`);
    return T;
}

/** Adjust values in a close-to-rotation matrix `a` to ensure it is a proper rotation matrix
 * (i.e. its columns and rows are orthonormal and determinant equal to 1, within available precission). */
function ensureRotationMatrix(out: Mat3, a: Mat3) {
    const x = Vec3.fromArray(_tmpVecX, a, 0);
    const y = Vec3.fromArray(_tmpVecY, a, 3);
    const z = Vec3.fromArray(_tmpVecZ, a, 6);
    Vec3.normalize(x, x);
    Vec3.orthogonalize(y, x, y);
    Vec3.normalize(z, Vec3.cross(z, x, y));
    Mat3.fromColumns(out, x, y, z);
    return out;
}
const _tmpVecX = Vec3();
const _tmpVecY = Vec3();
const _tmpVecZ = Vec3();

/** Create an array of props for `TransformStructureConformation` transformers from all 'transform' nodes applied to a 'structure' node. */
export function transformProps(node: SubTreeOfKind<MolstarTree, 'structure'>): StateTransformer.Params<TransformStructureConformation>[] {
    const result = [] as StateTransformer.Params<TransformStructureConformation>[];
    const transforms = getChildren(node).filter(c => c.kind === 'transform') as MolstarNode<'transform'>[];
    for (const transform of transforms) {
        const { rotation, translation } = transform.params;
        const matrix = transformFromRotationTranslation(rotation, translation);
        result.push({ transform: { name: 'matrix', params: { data: matrix, transpose: false } } });
    }
    return result;
}

/** Collect distinct annotation specs from all nodes in `tree` and set `context.annotationMap[node]` to respective annotationIds */
export function collectAnnotationReferences(tree: SubTree<MolstarTree>, context: MolstarLoadingContext): MVSAnnotationSpec[] {
    const distinctSpecs: { [key: string]: MVSAnnotationSpec } = {};
    dfs(tree, node => {
        let spec: Omit<MVSAnnotationSpec, 'id'> | undefined = undefined;
        if (AnnotationFromUriKinds.has(node.kind as any)) {
            const p = (node as MolstarNode<AnnotationFromUriKind>).params;
            spec = { source: { name: 'url', params: { url: p.uri, format: p.format } }, schema: p.schema, cifBlock: blockSpec(p.block_header, p.block_index), cifCategory: p.category_name ?? undefined };
        } else if (AnnotationFromSourceKinds.has(node.kind as any)) {
            const p = (node as MolstarNode<AnnotationFromSourceKind>).params;
            spec = { source: { name: 'source-cif', params: {} }, schema: p.schema, cifBlock: blockSpec(p.block_header, p.block_index), cifCategory: p.category_name ?? undefined };
        }
        if (spec) {
            const key = canonicalJsonString(spec as any);
            distinctSpecs[key] ??= { ...spec, id: stringHash(key) };
            (context.annotationMap ??= new Map()).set(node, distinctSpecs[key].id);
        }
    });
    return Object.values(distinctSpecs);
}
function blockSpec(header: string | null | undefined, index: number | null | undefined): MVSAnnotationSpec['cifBlock'] {
    if (isDefined(header)) {
        return { name: 'header', params: { header: header } };
    } else {
        return { name: 'index', params: { index: index ?? 0 } };
    }
}

/** Collect annotation tooltips from all nodes in `tree` and map them to annotationIds. */
export function collectAnnotationTooltips(tree: SubTreeOfKind<MolstarTree, 'structure'>, context: MolstarLoadingContext) {
    const annotationTooltips: MVSAnnotationTooltipsProps['tooltips'] = [];
    dfs(tree, node => {
        if (node.kind === 'tooltip_from_uri' || node.kind === 'tooltip_from_source') {
            const annotationId = context.annotationMap?.get(node);
            if (annotationId) {
                annotationTooltips.push({ annotationId, fieldName: node.params.field_name });
            };
        }
    });
    return arrayDistinct(annotationTooltips);
}
/** Collect annotation tooltips from all nodes in `tree`. */
export function collectInlineTooltips(tree: SubTreeOfKind<MolstarTree, 'structure'>, context: MolstarLoadingContext) {
    const inlineTooltips: CustomTooltipsProps['tooltips'] = [];
    dfs(tree, (node, parent) => {
        if (node.kind === 'tooltip') {
            if (parent?.kind === 'component') {
                inlineTooltips.push({
                    text: node.params.text,
                    selector: componentPropsFromSelector(parent.params.selector),
                });
            } else if (parent?.kind === 'component_from_uri' || parent?.kind === 'component_from_source') {
                const p = componentFromXProps(parent, context);
                if (isDefined(p.annotationId) && isDefined(p.fieldName) && isDefined(p.fieldValues)) {
                    inlineTooltips.push({
                        text: node.params.text,
                        selector: {
                            name: 'annotation',
                            params: { annotationId: p.annotationId, fieldName: p.fieldName, fieldValues: p.fieldValues },
                        },
                    });
                }
            }
        }
    });
    return inlineTooltips;
}

/** Return `true` for components nodes which only serve for tooltip placement (not to be created in the MolStar object hierarchy) */
export function isPhantomComponent(node: SubTreeOfKind<MolstarTree, 'component' | 'component_from_uri' | 'component_from_source'>) {
    return node.children && node.children.every(child => child.kind === 'tooltip' || child.kind === 'tooltip_from_uri' || child.kind === 'tooltip_from_source');
    // These nodes could theoretically be removed when converting MVS to Molstar tree, but would get very tricky if we allow nested components
}

/** Create props for `StructureFromModel` transformer from a structure node. */
export function structureProps(node: MolstarNode<'structure'>): StateTransformer.Params<StructureFromModel> {
    const params = node.params;
    switch (params.type) {
        case 'model':
            return {
                type: {
                    name: 'model',
                    params: {}
                },
            };
        case 'assembly':
            return {
                type: {
                    name: 'assembly',
                    params: { id: params.assembly_id ?? undefined }
                },
            };
        case 'symmetry':
            return {
                type: {
                    name: 'symmetry',
                    params: { ijkMin: Vec3.ofArray(params.ijk_min), ijkMax: Vec3.ofArray(params.ijk_max) }
                },
            };
        case 'symmetry_mates':
            return {
                type: {
                    name: 'symmetry-mates',
                    params: { radius: params.radius }
                }
            };
        default:
            throw new Error(`NotImplementedError: Loading action for "structure" node, type "${params.type}"`);
    }
}

/** Create value for `type` prop for `StructureComponent` transformer based on a MVS selector. */
export function componentPropsFromSelector(selector?: ParamsOfKind<MolstarTree, 'component'>['selector']): StructureComponentParams['type'] {
    if (selector === undefined) {
        return SelectorAll;
    } else if (typeof selector === 'string') {
        return { name: 'static', params: selector };
    } else if (Array.isArray(selector)) {
        return { name: 'expression', params: rowsToExpression(selector) };
    } else {
        return { name: 'expression', params: rowToExpression(selector) };
    }
}

/** Create props for `StructureRepresentation3D` transformer from a label_from_* node. */
export function labelFromXProps(node: MolstarNode<'label_from_uri' | 'label_from_source'>, context: MolstarLoadingContext): Partial<StateTransformer.Params<StructureRepresentation3D>> {
    const annotationId = context.annotationMap?.get(node);
    const fieldName = node.params.field_name;
    const nearestReprNode = context.nearestReprMap?.get(node);
    return {
        type: { name: MVSAnnotationLabelRepresentationProvider.name, params: { annotationId, fieldName } },
        colorTheme: colorThemeForNode(nearestReprNode, context),
    };
}

/** Create props for `AnnotationStructureComponent` transformer from a component_from_* node. */
export function componentFromXProps(node: MolstarNode<'component_from_uri' | 'component_from_source'>, context: MolstarLoadingContext): Partial<MVSAnnotationStructureComponentProps> {
    const annotationId = context.annotationMap?.get(node);
    const { field_name, field_values } = node.params;
    return {
        annotationId,
        fieldName: field_name,
        fieldValues: field_values ? { name: 'selected', params: field_values.map(v => ({ value: v })) } : { name: 'all', params: {} },
        nullIfEmpty: false,
    };
}

/** Create props for `StructureRepresentation3D` transformer from a representation node. */
export function representationProps(params: ParamsOfKind<MolstarTree, 'representation'>): Partial<StateTransformer.Params<StructureRepresentation3D>> {
    switch (params.type) {
        case 'cartoon':
            return {
                type: { name: 'cartoon', params: {} },
            };
        case 'ball_and_stick':
            return {
                type: { name: 'ball-and-stick', params: { sizeFactor: 0.5, sizeAspectRatio: 0.5 } },
            };
        case 'surface':
            return {
                type: { name: 'molecular-surface', params: {} },
                sizeTheme: { name: 'physical', params: { scale: 1 } },
            };
        default:
            throw new Error('NotImplementedError');
    }
}

/** Create value for `colorTheme` prop for `StructureRepresentation3D` transformer from a representation node based on color* nodes in its subtree. */
export function colorThemeForNode(node: SubTreeOfKind<MolstarTree, 'color' | 'color_from_uri' | 'color_from_source' | 'representation'> | undefined, context: MolstarLoadingContext): StateTransformer.Params<StructureRepresentation3D>['colorTheme'] {
    if (node?.kind === 'representation') {
        const children = getChildren(node).filter(c => c.kind === 'color' || c.kind === 'color_from_uri' || c.kind === 'color_from_source') as MolstarNode<'color' | 'color_from_uri' | 'color_from_source'>[];
        if (children.length === 0) {
            return {
                name: 'uniform',
                params: { value: decodeColor(DefaultColor) },
            };
        } else if (children.length === 1 && appliesColorToWholeRepr(children[0])) {
            return colorThemeForNode(children[0], context);
        } else {
            const layers: MultilayerColorThemeProps['layers'] = children.map(
                c => ({ theme: colorThemeForNode(c, context), selection: componentPropsFromSelector(c.kind === 'color' ? c.params.selector : undefined) })
            );
            return {
                name: MultilayerColorThemeName,
                params: { layers },
            };
        }
    }
    let annotationId: string | undefined = undefined;
    let fieldName: string | undefined = undefined;
    let color: string | undefined = undefined;
    switch (node?.kind) {
        case 'color_from_uri':
        case 'color_from_source':
            annotationId = context.annotationMap?.get(node);
            fieldName = node.params.field_name;
            break;
        case 'color':
            color = node.params.color;
            break;
    }
    if (annotationId) {
        return {
            name: MVSAnnotationColorThemeProvider.name,
            params: { annotationId, fieldName, background: NoColor } satisfies Partial<MVSAnnotationColorThemeProps>,
        };
    } else {
        return {
            name: 'uniform',
            params: { value: decodeColor(color) },
        };
    }
}
function appliesColorToWholeRepr(node: MolstarNode<'color' | 'color_from_uri' | 'color_from_source'>): boolean {
    if (node.kind === 'color') {
        return !isDefined(node.params.selector) || node.params.selector === 'all';
    } else {
        return true;
    }
}

/** Create a mapping of nearest representation nodes for each node in the tree
 * (to transfer coloring to label nodes smartly).
 * Only considers nodes within the same 'structure' subtree. */
export function makeNearestReprMap(root: MolstarTree) {
    const map = new Map<MolstarNode, MolstarNode<'representation'>>();
    // Propagate up:
    dfs(root, undefined, (node, parent) => {
        if (node.kind === 'representation') {
            map.set(node, node);
        }
        if (node.kind !== 'structure' && map.has(node) && parent && !map.has(parent)) { // do not propagate above the lowest structure node
            map.set(parent, map.get(node)!);
        }
    });
    // Propagate down:
    dfs(root, (node, parent) => {
        if (!map.has(node) && parent && map.has(parent)) {
            map.set(node, map.get(parent)!);
        }
    });
    return map;
}
