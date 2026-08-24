/**
 * Copyright (c) 2023-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Adam Midlik <midlik@gmail.com>
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Mat3, Mat4, Quat, Vec3 } from '../../mol-math/linear-algebra';
import { Volume } from '../../mol-model/volume';
import { StructureComponentParams } from '../../mol-plugin-state/helpers/structure-component';
import { StructureFromModel, StructureInstances, TransformStructureConformation } from '../../mol-plugin-state/transforms/model';
import { StructureRepresentation3D, VolumeRepresentation3D } from '../../mol-plugin-state/transforms/representation';
import { VolumeInstances, VolumeTransform } from '../../mol-plugin-state/transforms/volume';
import { StateTransformer } from '../../mol-state';
import { arrayDistinct } from '../../mol-util/array';
import { Clip } from '../../mol-util/clip';
import { canonicalJsonString } from '../../mol-util/json';
import { stringToWords } from '../../mol-util/string';
import { MVSAnnotationColorThemeProps, MVSAnnotationColorThemeProvider, MVSCategoricalPaletteProps, MVSContinuousPaletteProps, MVSDiscretePaletteProps } from './components/annotation-color-theme';
import { MVSAnnotationLabelProps, MVSAnnotationLabelRepresentationProvider } from './components/annotation-label/representation';
import { MVSAnnotationSpec } from './components/annotation-prop';
import { MVSAnnotationStructureComponentProps } from './components/annotation-structure-component';
import { MVSAnnotationTooltipsProps } from './components/annotation-tooltips-prop';
import { CustomLabelTextProps } from './components/custom-label/visual';
import { CustomTooltipsProps } from './components/custom-tooltips-prop';
import { MultilayerColorThemeName, MultilayerColorThemeProps, NoColor } from './components/multilayer-color-theme';
import { SelectorAll } from './components/selector';
import { MVSSplitUniformColorThemeProps, MVSSplitUniformColorThemeProvider, OptionalSplitColorProp, SplitColorProp } from './components/split-uniform-color-theme';
import { MvsNamedColorDicts, MvsNamedColorLists } from './helpers/colors';
import { rowToExpression, rowsToExpression } from './helpers/selections';
import { ElementOfSet, decodeColor, isDefined, stringHash } from './helpers/utils';
import { MolstarLoadingContext } from './load';
import { UpdateTarget, mvsRefTags } from './load-generic';
import { Subtree, getChildren } from './tree/generic/tree-schema';
import { dfs, formatObject } from './tree/generic/tree-utils';
import { MolstarKind, MolstarNode, MolstarNodeParams, MolstarSubtree, MolstarTree } from './tree/molstar/molstar-tree';
import { DefaultColor } from './tree/mvs/mvs-tree';
import { CategoricalPalette, CategoricalPaletteDefaults, ColorDictNameT, ColorListNameT, ColorT, ContinuousPalette, ContinuousPaletteDefaults, DiscretePalette, DiscretePaletteDefaults } from './tree/mvs/param-types';


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

export function decomposeRotationMatrix(rotation: number[] | null | undefined) {
    if (rotation && rotation.length !== 9) throw new Error(`'rotation' param for 'transform' node must be array of 9 elements, found ${rotation}`);
    if (rotation) {
        const rotMatrix = Mat3.fromArray(Mat3(), rotation, 0);
        ensureRotationMatrix(rotMatrix, rotMatrix);
        const quat = Quat.fromMat3(Quat(), rotMatrix);
        const axis = Vec3();
        const angle = Quat.getAxisAngle(axis, quat) * 180 / Math.PI;
        return { axis, angle };
    }
    return { axis: Vec3.create(1, 0, 0), angle: 0 };
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

export function transformAndInstantiateStructure(
    target: UpdateTarget,
    node: MolstarSubtree<'structure' | 'component' | 'component_from_source' | 'component_from_uri'>,
) {
    return applyTransformAndInstances(target, node, TransformStructureConformation, StructureInstances);
}

export function transformAndInstantiateVolume(target: UpdateTarget, node: MolstarSubtree<'volume'>) {
    return applyTransformAndInstances(target, node, VolumeTransform, VolumeInstances);
}

function applyTransformAndInstances(target: UpdateTarget, node: MolstarSubtree, transform: StateTransformer, instantiate: StateTransformer) {
    let modified = target;
    for (const { params, ref } of transformProps(node, 'transform')) {
        modified = UpdateTarget.apply(modified, transform, params);
        UpdateTarget.tag(modified, mvsRefTags(ref));
    }

    const instances = transformProps(node, 'instance');
    if (instances.length > 0) {
        modified = UpdateTarget.apply(modified, instantiate, { transforms: instances.map(i => i.params) });
    }

    return modified;
}

/** Create an array of props for `TransformStructureConformation` transformers from all 'transform' nodes applied to a 'structure' node. */
function transformProps(node: MolstarSubtree, kind: 'transform' | 'instance') {
    const result = [] as { params: StateTransformer.Params<TransformStructureConformation>, ref?: string }[];
    const transforms = getChildren(node).filter(c => c.kind === kind) as MolstarNode<'transform'>[];
    for (const transform of transforms) {
        let matrix: Mat4 | undefined = transform.params.matrix as Mat4 | undefined;
        if (!matrix) {
            const { rotation, translation, rotation_center } = transform.params;
            if (rotation_center) {
                const axisAngle = decomposeRotationMatrix(rotation);
                result.push({
                    params: {
                        transform: {
                            name: 'components',
                            params: {
                                translation: translation ? Vec3.fromArray(Vec3(), translation, 0) : Vec3.create(0, 0, 0),
                                angle: axisAngle.angle,
                                axis: axisAngle.axis,
                                rotationCenter: rotation_center === 'centroid'
                                    ? { name: 'centroid', params: {} }
                                    : { name: 'point', params: { point: Vec3.fromArray(Vec3(), rotation_center, 0) } }
                            }
                        }
                    },
                    ref: transform.ref
                });
                continue;
            }
            matrix = transformFromRotationTranslation(rotation, translation);
        }
        result.push({ params: { transform: { name: 'matrix', params: { data: matrix, transpose: false } } }, ref: transform.ref });
    }

    return result;
}

/** Collect distinct annotation specs from all nodes in `tree` and set `context.annotationMap[node]` to respective annotationIds */
export function collectAnnotationReferences(tree: Subtree<MolstarTree>, context: MolstarLoadingContext): MVSAnnotationSpec[] {
    const distinctSpecs: { [key: string]: MVSAnnotationSpec } = {};
    dfs(tree, node => {
        let spec: Omit<MVSAnnotationSpec, 'id'> | undefined = undefined;
        if (AnnotationFromUriKinds.has(node.kind as any)) {
            const p = (node as MolstarNode<AnnotationFromUriKind>).params;
            spec = {
                source: { name: 'url', params: { url: p.uri, format: p.format } }, schema: p.schema,
                cifBlock: blockSpec(p.block_header, p.block_index), cifCategory: p.category_name,
                fieldRemapping: Object.entries(p.field_remapping).map(([key, value]) => ({ standardName: key, actualName: value })),
            };
        } else if (AnnotationFromSourceKinds.has(node.kind as any)) {
            const p = (node as MolstarNode<AnnotationFromSourceKind>).params;
            spec = {
                source: { name: 'source-cif', params: {} }, schema: p.schema,
                cifBlock: blockSpec(p.block_header, p.block_index), cifCategory: p.category_name,
                fieldRemapping: Object.entries(p.field_remapping).map(([key, value]) => ({ standardName: key, actualName: value })),
            };
        }
        if (spec) {
            const key = canonicalJsonString(spec as any);
            distinctSpecs[key] ??= { ...spec, id: stringHash(key) };
            context.annotationMap.set(node as MolstarNode<AnnotationFromUriKind | AnnotationFromSourceKind>, distinctSpecs[key].id);
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
export function collectAnnotationTooltips(tree: MolstarSubtree<'structure'>, context: MolstarLoadingContext): MVSAnnotationTooltipsProps['tooltips'] {
    const annotationTooltips: MVSAnnotationTooltipsProps['tooltips'] = [];
    dfs(tree, node => {
        if (node.kind === 'tooltip_from_uri' || node.kind === 'tooltip_from_source') {
            const annotationId = context.annotationMap.get(node);
            if (annotationId) {
                annotationTooltips.push({ annotationId, fieldName: node.params.field_name, textFormat: node.params.text_format });
            };
        }
    });
    return arrayDistinct(annotationTooltips);
}
/** Collect inline tooltips from all nodes in `tree`. */
export function collectInlineTooltips(tree: MolstarSubtree<'structure'>, context: MolstarLoadingContext): CustomTooltipsProps['tooltips'] {
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
                            params: { annotationId: p.annotationId, fieldName: p.fieldName, fieldValues: p.fieldValues, label: p.label || 'Annotation' },
                        },
                    });
                }
            }
        }
    });
    return inlineTooltips;
}
/** Collect inline labels from all nodes in `tree`. */
export function collectInlineLabels(tree: MolstarSubtree<'structure'>, context: MolstarLoadingContext): CustomLabelTextProps['items'] {
    const inlineLabels: CustomLabelTextProps['items'] = [];
    dfs(tree, (node, parent) => {
        if (node.kind === 'label') {
            if (parent?.kind === 'component') {
                inlineLabels.push({
                    text: node.params.text,
                    position: {
                        name: 'selection',
                        params: {
                            selector: componentPropsFromSelector(parent.params.selector),
                        },
                    },
                });
            } else if (parent?.kind === 'component_from_uri' || parent?.kind === 'component_from_source') {
                const p = componentFromXProps(parent, context);
                if (isDefined(p.annotationId) && isDefined(p.fieldName) && isDefined(p.fieldValues)) {
                    inlineLabels.push({
                        text: node.params.text,
                        position: {
                            name: 'selection',
                            params: {
                                selector: {
                                    name: 'annotation',
                                    params: { annotationId: p.annotationId, fieldName: p.fieldName, fieldValues: p.fieldValues, label: p.label || 'Annotation' },
                                },
                            },
                        },
                    });
                }
            }
        }
    });
    return inlineLabels;
}

/** Return `true` for components nodes which only serve for tooltip placement (not to be created in the MolStar object hierarchy) */
export function isPhantomComponent(node: MolstarSubtree<'component' | 'component_from_uri' | 'component_from_source'>) {
    if (node.ref !== undefined) return false;
    if (node.custom !== undefined && Object.keys(node.custom).length > 0) return false;
    return node.children && node.children.every(child => child.kind === 'tooltip' || child.kind === 'label');
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
export function componentPropsFromSelector(selector?: MolstarNodeParams<'component'>['selector']): StructureComponentParams['type'] {
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

/** Return a pretty name for a value of selector param, e.g.  "protein" -> 'Protein', {label_asym_id: "A"} -> 'Custom Selection: {label_asym_id: "A"}' */
export function prettyNameFromSelector(selector?: MolstarNodeParams<'component'>['selector']): string {
    if (selector === undefined) {
        return 'All';
    } else if (typeof selector === 'string') {
        return stringToWords(selector);
    } else if (Array.isArray(selector)) {
        return `Custom Selection: [${selector.map(formatObject).join(', ')}]`;
    } else {
        return `Custom Selection: ${formatObject(selector)}`;
    }
}

/** Create props for `StructureRepresentation3D` transformer from a label_from_* node. */
export function labelFromXProps(node: MolstarNode<'label_from_uri' | 'label_from_source'>, context: MolstarLoadingContext): Partial<StateTransformer.Params<StructureRepresentation3D>> {
    const annotationId = context.annotationMap.get(node)!;
    const fieldName = node.params.field_name;
    const textFormat = node.params.text_format;
    const groupBy = node.params.group_by_fields ?? [node.params.field_remapping.group_id ?? 'group_id'];
    const nearestReprNode = context.nearestReprMap?.get(node);
    return {
        type: { name: MVSAnnotationLabelRepresentationProvider.name, params: { annotationId, fieldName, textFormat, groupByFields: groupBy.map(x => ({ fieldName: x })) } satisfies Partial<MVSAnnotationLabelProps> },
        colorTheme: colorThemeForNode(nearestReprNode, context),
    };
}

/** Create props for `AnnotationStructureComponent` transformer from a component_from_* node. */
export function componentFromXProps(node: MolstarNode<'component_from_uri' | 'component_from_source'>, context: MolstarLoadingContext): Partial<MVSAnnotationStructureComponentProps> {
    const annotationId = context.annotationMap.get(node);
    const { field_name, field_values } = node.params;
    return {
        annotationId,
        fieldName: field_name,
        fieldValues: field_values ? { name: 'selected', params: field_values.map(v => ({ value: v })) } : { name: 'all', params: {} },
        nullIfEmpty: false,
    };
}

/** Create props for `StructureRepresentation3D` transformer from a representation node. */
function representationPropsBase(node: MolstarSubtree<'representation'>): Partial<StateTransformer.Params<StructureRepresentation3D>> {
    const alpha = alphaForNode(node);
    const params = node.params;
    switch (params.type) {
        case 'cartoon':
            return {
                type: { name: 'cartoon', params: { alpha, tubularHelices: params.tubular_helices } },
                sizeTheme: { name: 'uniform', params: { value: params.size_factor } },
            };
        case 'backbone':
            return {
                type: { name: 'backbone', params: { alpha } },
                sizeTheme: { name: 'uniform', params: { value: params.size_factor } },
            };
        case 'ball_and_stick':
            return {
                type: { name: 'ball-and-stick', params: { sizeFactor: 0.5, sizeAspectRatio: 0.5, alpha, ignoreHydrogens: params.ignore_hydrogens } },
                sizeTheme: { name: 'uniform', params: { value: params.size_factor } },
            };
        case 'line':
            return {
                type: { name: 'line', params: { alpha, ignoreHydrogens: params.ignore_hydrogens } },
                sizeTheme: { name: 'uniform', params: { value: params.size_factor } },
            };
        case 'spacefill':
            return {
                type: { name: 'spacefill', params: { alpha, ignoreHydrogens: params.ignore_hydrogens } },
                sizeTheme: { name: 'physical', params: { scale: params.size_factor } },
            };
        case 'carbohydrate':
            return {
                type: { name: 'carbohydrate', params: { alpha, sizeFactor: 1.75 } },
                sizeTheme: { name: 'uniform', params: { value: params.size_factor } },
            };
        case 'surface': {
            return {
                type: {
                    name: params.surface_type === 'gaussian' ? 'gaussian-surface' : 'molecular-surface',
                    params: { alpha, ignoreHydrogens: params.ignore_hydrogens }
                },
                sizeTheme: { name: 'physical', params: { scale: params.size_factor } },
            };
        }
        case 'putty': {
            const sizeTheme = params.size_theme ?? 'uniform';
            return {
                type: { name: 'putty', params: { alpha, sizeFactor: params.size_factor } },
                sizeTheme: sizeTheme === 'uncertainty'
                    ? { name: 'uncertainty', params: {} }
                    : { name: 'uniform', params: { value: params.size_factor } },
            };
        }
        default:
            throw new Error('NotImplementedError');
    }
}

export function representationProps(node: MolstarSubtree<'representation'>): Partial<StateTransformer.Params<StructureRepresentation3D>> {
    const base = representationPropsBase(node);
    const clip = clippingForNode(node);
    if (clip) {
        base.type!.params = { ...base.type?.params, clip };
    }
    if (node.custom?.molstar_representation_params) {
        base.type!.params = { ...base.type!.params, ...node.custom.molstar_representation_params };
    }
    return base;
}

/** Create value for `type.params.alpha` prop for `StructureRepresentation3D` transformer from a representation node based on 'opacity' nodes in its subtree. */
export function alphaForNode(node: MolstarSubtree<'representation' | 'volume_representation' | 'shape'>): number {
    const children = getChildren(node).filter(c => c.kind === 'opacity');
    if (children.length > 0) {
        return children[children.length - 1].params.opacity;
    } else {
        return 1;
    }
}

function getCommonClipParams(node: MolstarNode<'clip'>): Pick<Clip.Props['objects'][number], 'invert' | 'transform'> {
    return {
        invert: !!node.params.invert,
        transform: node.params.check_transform ? Mat4.fromArray(Mat4(), node.params.check_transform, 0) : Mat4.identity(),
    };
}

function getClipObject(node: MolstarNode<'clip'>): Clip.Props['objects'][number] | undefined {
    switch (node.params.type) {
        case 'sphere':
            return {
                type: 'sphere',
                position: Vec3.ofArray(node.params.center),
                scale: typeof node.params.radius === 'number'
                    ? Vec3.create(2 * node.params.radius, 2 * node.params.radius, 2 * node.params.radius)
                    : Vec3.create(2, 2, 2),
                rotation: { axis: Vec3.create(1, 0, 0), angle: 0 },
                ...getCommonClipParams(node),
            };
        case 'plane': {
            const up = Vec3.create(0, 1, 0);
            const n = Vec3.normalize(Vec3(), Vec3.ofArray(node.params.normal));
            const axis = Vec3.cross(Vec3(), up, n);
            const isSingular = Vec3.magnitude(axis) < 1e-6;
            return {
                type: 'plane',
                position: Vec3.ofArray(node.params.point),
                scale: Vec3.create(1, 1, 1),
                rotation: {
                    axis: isSingular ? Vec3.unitX : axis,
                    angle: isSingular ? 0 : Vec3.angle(up, n) * 180 / Math.PI,
                },
                ...getCommonClipParams(node),
            };
        }
        case 'box':
            const q = Quat.fromMat3(Quat(), Mat3.fromArray(Mat3(), node.params.rotation, 0));
            const axis = Vec3();
            const angle = Quat.getAxisAngle(axis, q) * 180 / Math.PI;
            return {
                type: 'cube',
                position: Vec3.ofArray(node.params.center),
                scale: Vec3.ofArray(node.params.size),
                rotation: { axis, angle },
                ...getCommonClipParams(node),
            };
        default:
            console.warn(`Mol* MVS: Unsupported clip type "${(node as MolstarNode<'clip'>).params.type}" in node ${node.ref}.`);
    }
}

export function clippingForNode(node: MolstarSubtree<'representation' | 'volume_representation' | 'primitives' | 'primitives_from_uri' | 'shape'>): Clip.Props | undefined {
    const children = getChildren(node).filter(c => c.kind === 'clip');
    if (!children.length) return;

    const variant = children[0].params.variant === 'object' ? 'instance' : 'pixel';
    const objects: Clip.Props['objects'] = children.map(getClipObject).filter(o => !!o);

    return { variant, objects } satisfies Clip.Props;
}

function hasMolStarUseDefaultColoring(node: MolstarNode): boolean {
    if (!node.custom) return false;
    return 'molstar_use_default_coloring' in node.custom || 'molstar_color_theme_name' in node.custom;
}

function customColoring(custom: any) {
    if (custom?.molstar_use_default_coloring) return undefined;
    return {
        name: custom?.molstar_color_theme_name ?? undefined,
        params: custom?.molstar_color_theme_params ?? {},
    };
}

/** Create value for `colorTheme` prop for `StructureRepresentation3D` transformer from a representation node based on color* nodes in its subtree. */
export function colorThemeForNode(node: MolstarSubtree<'color' | 'color_from_uri' | 'color_from_source' | 'representation' | 'volume'> | undefined, context: MolstarLoadingContext): StateTransformer.Params<StructureRepresentation3D>['colorTheme'] | undefined {
    if (node?.kind === 'representation') {
        const children = getChildren(node).filter(c => c.kind === 'color' || c.kind === 'color_from_uri' || c.kind === 'color_from_source') as MolstarNode<'color' | 'color_from_uri' | 'color_from_source'>[];
        if (children.length === 0) {
            return uniformColorTheme(DefaultColor);
        } else if (children.length === 1 && hasMolStarUseDefaultColoring(children[0])) {
            return customColoring(children[0].custom);
        } else if (children.length === 1 && appliesColorToWholeRepr(children[0])) {
            return colorThemeForNode(children[0], context);
        } else {
            const layers: MultilayerColorThemeProps['layers'] = children.map(
                c => {
                    const theme = colorThemeForNode(c, context);
                    if (!theme) return undefined;
                    return { theme, selection: componentPropsFromSelector(c.params.selector) };
                }
            ).filter(t => !!t);
            return {
                name: MultilayerColorThemeName,
                params: { layers },
            };
        }
    }
    if (node?.kind === 'color') {
        if (hasMolStarUseDefaultColoring(node)) {
            return customColoring(node.custom);
        }
        return uniformColorTheme(node.params.color);
    }
    if (node?.kind === 'color_from_uri' || node?.kind === 'color_from_source') {
        const annotationId = context.annotationMap.get(node);
        if (annotationId === undefined) return {
            name: 'uniform',
            params: {},
        };

        const fieldName = node.params.field_name;
        return {
            name: MVSAnnotationColorThemeProvider.name,
            params: { annotationId, fieldName, background: NoColor, palette: palettePropsFromMVSPalette(node.params.palette) } satisfies Partial<MVSAnnotationColorThemeProps>,
        };
    }
}

function appliesColorToWholeRepr(node: MolstarNode<'color' | 'color_from_uri' | 'color_from_source'>): boolean {
    return !isDefined(node.params.selector) || node.params.selector === 'all';
}

/** Create parameters for `uniform` (or `mvs-split-uniform`) color theme, supports split coloring. */
function uniformColorTheme(color: ColorT): StateTransformer.Params<StructureRepresentation3D>['colorTheme'] {
    if (color.includes('/')) {
        return {
            name: MVSSplitUniformColorThemeProvider.name,
            params: { value: SplitColorProp.fromString(color) } satisfies MVSSplitUniformColorThemeProps,
        };
    } else {
        return {
            name: 'uniform',
            params: { value: decodeColor(color) },
        };
    }
}

export function palettePropsFromMVSPalette(palette: MolstarNode<'color_from_uri' | 'color_from_source'>['params']['palette']): MVSAnnotationColorThemeProps['palette'] {
    if (!palette) {
        return { name: 'direct', params: {} };
    }
    if (palette.kind === 'categorical') {
        const fullParams: Required<CategoricalPalette> = objMerge(CategoricalPaletteDefaults, palette);
        return {
            name: 'categorical',
            params: {
                colors: categoricalPalettePropsFromMVSColors(fullParams.colors),
                repeatColorList: fullParams.repeat_color_list,
                sort: fullParams.sort,
                sortDirection: fullParams.sort_direction,
                caseInsensitive: fullParams.case_insensitive,
                missingColor: SplitColorProp.fromString(fullParams.missing_color),
            } satisfies MVSCategoricalPaletteProps,
        };
    }
    if (palette.kind === 'discrete') {
        const fullParams: Required<DiscretePalette> = objMerge(DiscretePaletteDefaults, palette);
        return {
            name: 'discrete',
            params: {
                colors: discretePalettePropsFromMVSColors(fullParams.colors, fullParams.reverse),
                mode: fullParams.mode,
                xMin: fullParams.value_domain[0],
                xMax: fullParams.value_domain[1],
            } satisfies MVSDiscretePaletteProps,
        };
    }
    if (palette.kind === 'continuous') {
        const fullParams: Required<ContinuousPalette> = objMerge(ContinuousPaletteDefaults, palette);
        const colors = continuousPalettePropsFromMVSColors(fullParams.colors, fullParams.reverse);
        return {
            name: 'continuous',
            params: {
                colors: colors,
                mode: fullParams.mode,
                xMin: fullParams.value_domain[0],
                xMax: fullParams.value_domain[1],
                underflowColor: fullParams.underflow_color === 'auto' ? minColor(colors) : SplitColorProp.fromString(fullParams.underflow_color),
                overflowColor: fullParams.overflow_color === 'auto' ? maxColor(colors) : SplitColorProp.fromString(fullParams.overflow_color),
            } satisfies MVSContinuousPaletteProps,
        };
    }
    throw new Error(`NotImplementedError: palettePropsFromMVSPalette is not implemented for palette kind "${(palette as any).kind}"`);
}

/** Merge properties of two object into a new object. Property values from `second` override those from `first`, but `undefined` is treated as if property missing while `null` as a regular value. */
function objMerge<T extends object, U extends object>(first: T, second: U): T & U {
    const out: Partial<T & U> = { ...first };
    for (const key in second) {
        const value = second[key];
        if (value !== undefined) out[key] = value as any;
    }
    return out as T & U;
}

function categoricalPalettePropsFromMVSColors(colors: Required<CategoricalPalette>['colors']): MVSCategoricalPaletteProps['colors'] {
    if (typeof colors === 'string') {
        // Named color list or dict
        if (colors in MvsNamedColorLists) {
            colors = MvsNamedColorLists[colors as ColorListNameT];
        } else if (colors in MvsNamedColorDicts) {
            colors = MvsNamedColorDicts[colors as ColorDictNameT];
        } else {
            console.warn(`Could not find named color palette "${colors}"`);
            return { name: 'list', params: [] };
        }
    }

    if (Array.isArray(colors)) {
        // Color list
        return { name: 'list', params: colors.map(c => ({ color: SplitColorProp.fromString(c) })) };
    } else {
        // Color dict
        return { name: 'dictionary', params: Object.entries(colors).map(([value, color]) => ({ value, color: SplitColorProp.fromString(color) })) };
    }
}


function discretePalettePropsFromMVSColors(colors: Required<DiscretePalette>['colors'], reverse: boolean): MVSDiscretePaletteProps['colors'] {
    if (typeof colors === 'string') {
        // Named color list
        if (colors in MvsNamedColorLists) {
            colors = MvsNamedColorLists[colors as ColorListNameT];
        } else {
            console.warn(`Could not find named color palette "${colors}"`);
            return [];
        }
    }

    if (colors.every(t => typeof t === 'string')) {
        // Color list without checkpoints
        if (reverse) colors = colors.slice().reverse();
        const sectionLength = 1 / colors.length;
        return colors.map((col, i) => ({ color: SplitColorProp.fromString(col), fromValue: i * sectionLength, toValue: (i + 1) * sectionLength }));
    } else if (colors.every(t => Array.isArray(t) && t.length === 2)) {
        // Color list with 1 checkpoint
        return colors.map(([col, from], i, arr) => ({ color: SplitColorProp.fromString(col), fromValue: from, toValue: arr[i + 1]?.[1] ?? Infinity }));
    } else if (colors.every(t => Array.isArray(t) && t.length === 3)) {
        // Color list with 2 checkpoint
        return colors.map(([col, from, to]) => ({ color: SplitColorProp.fromString(col), fromValue: from ?? -Infinity, toValue: to ?? Infinity }));
    } else {
        return []; // should be unreachable
    }
}

function continuousPalettePropsFromMVSColors(colors: Required<ContinuousPalette>['colors'], reverse: boolean): MVSContinuousPaletteProps['colors'] {
    if (typeof colors === 'string') {
        // Named color list
        if (colors in MvsNamedColorLists) {
            colors = MvsNamedColorLists[colors as ColorListNameT];
        } else {
            console.warn(`Could not find named color palette "${colors}"`);
            return [];
        }
    }

    if (colors.every(t => typeof t === 'string')) {
        // Color list without checkpoints
        if (reverse) colors = colors.slice().reverse();
        const n = colors.length - 1;
        return colors.map((col, i) => ({ color: SplitColorProp.fromString(col), checkpoint: i / n }));
    } else {
        // Color list with checkpoints
        // Not applying `reverse` here, as it would have no effect
        return colors.map(([col, checkpoint]) => ({ color: SplitColorProp.fromString(col), checkpoint: checkpoint }));
    }
}

/** Return the color with the lowest checkpoint. */
function minColor(colors: MVSContinuousPaletteProps['colors']): OptionalSplitColorProp {
    if (colors.length === 0) return SplitColorProp.none();
    return colors.reduce((a, b) => a.checkpoint < b.checkpoint ? a : b).color;
}
/** Return the color with the highest checkpoint. */
function maxColor(colors: MVSContinuousPaletteProps['colors']): OptionalSplitColorProp {
    if (colors.length === 0) return SplitColorProp.none();
    return colors.reduce((a, b) => a.checkpoint > b.checkpoint ? a : b).color;
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

/** Create props for `VolumeRepresentation3D` transformer from a representation node. */
export function volumeRepresentationProps(node: MolstarSubtree<'volume_representation'>): Partial<StateTransformer.Params<VolumeRepresentation3D>> {
    const alpha = alphaForNode(node);
    const clip = clippingForNode(node);
    const params = node.params;

    const isoValue = typeof params.absolute_isovalue === 'number' ? Volume.IsoValue.absolute(params.absolute_isovalue) : Volume.IsoValue.relative(params.relative_isovalue ?? 0);
    switch (params.type) {
        case 'isosurface':
            const visuals: ('wireframe' | 'solid')[] = [];
            if (params.show_wireframe) visuals.push('wireframe');
            if (params.show_faces) visuals.push('solid');
            return {
                type: { name: 'isosurface', params: { alpha, isoValue, visuals, clip } },
            };
        case 'grid_slice':
            const isRelative = params.relative_index !== undefined;
            const dimension = {
                name: isRelative ? `relative${params.dimension.toUpperCase()}` : params.dimension,
                params: params.relative_index ?? params.relative_index
            };
            return {
                type: { name: 'slice', params: { alpha, dimension, isoValue, clip } },
            };
        default:
            throw new Error('NotImplementedError');
    }
}

/** Create value for `colorTheme` prop for `StructureRepresentation3D` transformer from a representation node based on color* nodes in its subtree. */
export function volumeColorThemeForNode(node: MolstarSubtree<'volume_representation'> | undefined, context: MolstarLoadingContext): StateTransformer.Params<VolumeRepresentation3D>['colorTheme'] | undefined {
    if (node?.kind !== 'volume_representation') return undefined;

    const children = getChildren(node).filter(c => c.kind === 'color') as MolstarNode<'color'>[];
    if (children.length === 0) {
        return {
            name: 'uniform',
            params: { value: decodeColor(DefaultColor) },
        };
    } if (children.length === 1) {
        return colorThemeForNode(children[0], context);
    }
}
