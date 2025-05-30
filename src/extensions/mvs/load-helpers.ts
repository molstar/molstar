/**
 * Copyright (c) 2023-2025 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Adam Midlik <midlik@gmail.com>
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Mat3, Mat4, Vec3 } from '../../mol-math/linear-algebra';
import { Volume } from '../../mol-model/volume';
import { StructureComponentParams } from '../../mol-plugin-state/helpers/structure-component';
import { StructureFromModel, TransformStructureConformation } from '../../mol-plugin-state/transforms/model';
import { StructureRepresentation3D, VolumeRepresentation3D } from '../../mol-plugin-state/transforms/representation';
import { StateTransformer } from '../../mol-state';
import { ElementSymbolColors } from '../../mol-theme/color/element-symbol';
import { ResidueNameColors } from '../../mol-theme/color/residue-name';
import { arrayDistinct } from '../../mol-util/array';
import { Color } from '../../mol-util/color';
import { ColorListName, ColorLists } from '../../mol-util/color/lists';
import { canonicalJsonString } from '../../mol-util/json';
import { omitObjectKeys } from '../../mol-util/object';
import { stringToWords } from '../../mol-util/string';
import { MVSAnnotationColorThemeProps, MVSAnnotationColorThemeProvider, MVSCategoricalPaletteProps } from './components/annotation-color-theme';
import { MVSAnnotationLabelRepresentationProvider } from './components/annotation-label/representation';
import { MVSAnnotationSpec } from './components/annotation-prop';
import { MVSAnnotationStructureComponentProps } from './components/annotation-structure-component';
import { MVSAnnotationTooltipsProps } from './components/annotation-tooltips-prop';
import { CustomLabelTextProps } from './components/custom-label/visual';
import { CustomTooltipsProps } from './components/custom-tooltips-prop';
import { MultilayerColorThemeName, MultilayerColorThemeProps, NoColor } from './components/multilayer-color-theme';
import { SelectorAll } from './components/selector';
import { rowToExpression, rowsToExpression } from './helpers/selections';
import { ElementOfSet, decodeColor, isDefined, stringHash } from './helpers/utils';
import { MolstarLoadingContext } from './load';
import { Subtree, getChildren } from './tree/generic/tree-schema';
import { dfs, formatObject } from './tree/generic/tree-utils';
import { MolstarKind, MolstarNode, MolstarNodeParams, MolstarSubtree, MolstarTree } from './tree/molstar/molstar-tree';
import { DefaultColor } from './tree/mvs/mvs-tree';
import { CategoricalPalette, ColorDictNameT, ColorListNameT } from './tree/mvs/param-types';


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
export function transformProps(node: MolstarSubtree<'structure'>): StateTransformer.Params<TransformStructureConformation>[] {
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
export function collectAnnotationReferences(tree: Subtree<MolstarTree>, context: MolstarLoadingContext): MVSAnnotationSpec[] {
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
                annotationTooltips.push({ annotationId, fieldName: node.params.field_name });
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
                            params: { annotationId: p.annotationId, fieldName: p.fieldName, fieldValues: p.fieldValues },
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
                                    params: { annotationId: p.annotationId, fieldName: p.fieldName, fieldValues: p.fieldValues },
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
    const annotationId = context.annotationMap.get(node);
    const fieldName = node.params.field_name;
    const nearestReprNode = context.nearestReprMap?.get(node);
    return {
        type: { name: MVSAnnotationLabelRepresentationProvider.name, params: { annotationId, fieldName } },
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
export function representationProps(node: MolstarSubtree<'representation'>): Partial<StateTransformer.Params<StructureRepresentation3D>> {
    const alpha = alphaForNode(node);
    const params = node.params;
    switch (params.type) {
        case 'cartoon':
            return {
                type: { name: 'cartoon', params: { alpha, tubularHelices: params.tubular_helices } },
                sizeTheme: { name: 'uniform', params: { value: params.size_factor } },
            };
        case 'ball_and_stick':
            return {
                type: { name: 'ball-and-stick', params: { sizeFactor: (params.size_factor ?? 1) * 0.5, sizeAspectRatio: 0.5, alpha, ignoreHydrogens: params.ignore_hydrogens } },
            };
        case 'spacefill':
            return {
                type: { name: 'spacefill', params: { alpha, ignoreHydrogens: params.ignore_hydrogens } },
                sizeTheme: { name: 'physical', params: { scale: params.size_factor } },
            };
        case 'carbohydrate':
            return {
                type: { name: 'carbohydrate', params: { alpha, sizeFactor: params.size_factor ?? 1 } },
            };
        case 'surface':
            return {
                type: { name: 'molecular-surface', params: { alpha, ignoreHydrogens: params.ignore_hydrogens } },
                sizeTheme: { name: 'physical', params: { scale: params.size_factor } },
            };
        default:
            throw new Error('NotImplementedError');
    }
}

/** Create value for `type.params.alpha` prop for `StructureRepresentation3D` transformer from a representation node based on 'opacity' nodes in its subtree. */
export function alphaForNode(node: MolstarSubtree<'representation' | 'volume_representation'>): number {
    const children = getChildren(node).filter(c => c.kind === 'opacity');
    if (children.length > 0) {
        return children[children.length - 1].params.opacity;
    } else {
        return 1;
    }
}

function hasMolStarUseDefaultColoring(node: MolstarNode): boolean {
    if (!node.custom) return false;
    return 'molstar_use_default_coloring' in node.custom || 'molstar_color_theme_name' in node.custom;
}

/** Create value for `colorTheme` prop for `StructureRepresentation3D` transformer from a representation node based on color* nodes in its subtree. */
export function colorThemeForNode(node: MolstarSubtree<'color' | 'color_from_uri' | 'color_from_source' | 'representation'> | undefined, context: MolstarLoadingContext): StateTransformer.Params<StructureRepresentation3D>['colorTheme'] | undefined {
    if (node?.kind === 'representation') {
        const children = getChildren(node).filter(c => c.kind === 'color' || c.kind === 'color_from_uri' || c.kind === 'color_from_source') as MolstarNode<'color' | 'color_from_uri' | 'color_from_source'>[];
        if (children.length === 0) {
            return {
                name: 'uniform',
                params: { value: decodeColor(DefaultColor) },
            };
        } else if (children.length === 1 && hasMolStarUseDefaultColoring(children[0])) {
            if (children[0].custom?.molstar_use_default_coloring) return undefined;
            const custom = children[0].custom;
            return {
                name: custom?.molstar_color_theme_name ?? undefined,
                params: custom?.molstar_color_theme_params ?? {},
            };
        } else if (children.length === 1 && appliesColorToWholeRepr(children[0])) {
            return colorThemeForNode(children[0], context);
        } else {
            const layers: MultilayerColorThemeProps['layers'] = children.map(
                c => {
                    const theme = colorThemeForNode(c, context);
                    if (!theme) return undefined;
                    return { theme, selection: componentPropsFromSelector(c.kind === 'color' ? c.params.selector : undefined) };
                }
            ).filter(t => !!t);
            return {
                name: MultilayerColorThemeName,
                params: { layers },
            };
        }
    }
    if (node?.kind === 'color') {
        return {
            name: 'uniform',
            params: { value: decodeColor(node.params.color) },
        };
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
    if (node.kind === 'color') {
        return !isDefined(node.params.selector) || node.params.selector === 'all';
    } else {
        return true;
    }
}

const FALLBACK_COLOR = decodeColor(DefaultColor)!;

function palettePropsFromMVSPalette(palette: MolstarNode<'color_from_uri' | 'color_from_source'>['params']['palette']): MVSAnnotationColorThemeProps['palette'] {
    if (!palette) {
        return { name: 'direct', params: {} };
    }
    if (palette.kind === 'categorical') {
        return {
            name: 'categorical',
            params: {
                colors: categoricalPalettePropsFromMVSColors(palette.colors ?? 'Category10'),
                repeatColorList: palette.repeat_color_list ?? false,
                sort: palette.sort ?? 'none',
                sortDirection: palette.sort_direction ?? 'ascending',
                setMissingColor: !!palette.missing_color,
                missingColor: decodeColor(palette.missing_color) ?? FALLBACK_COLOR,
                // TODO specify defaults with param-types or mvs-tree
            },
        };
    }
    throw new Error(`NotImplementedError: palettePropsFromMVSPalette is not implemented for palette kind "${palette.kind}"`);
}

function categoricalPalettePropsFromMVSColors(colors: CategoricalPalette['colors']): MVSCategoricalPaletteProps['colors'] {
    if (typeof colors === 'string') {
        if (colors in MvsNamedColorListToMolstarName) {
            const molstarColorListName = MvsNamedColorListToMolstarName[colors as keyof typeof MvsNamedColorListToMolstarName];
            const colorList = ColorLists[molstarColorListName];
            if (colorList) {
                return { name: 'list', params: { kind: 'set', colors: colorList.list } };
            }
        }
        if (colors in MvsNamedColorDicts) {
            const dict = MvsNamedColorDicts[colors as keyof typeof MvsNamedColorDicts];
            return { name: 'dictionary', params: Object.entries(dict).map(([value, color]) => ({ value, color })) };
        }
        console.warn(`Could not find named color palette "${colors}"`);
    }
    if (Array.isArray(colors)) {
        return { name: 'list', params: { kind: 'set', colors: colors.map(c => decodeColor(c) ?? FALLBACK_COLOR) } };
    }
    if (typeof colors === 'object') {
        return { name: 'dictionary', params: Object.entries(colors).map(([value, color]) => ({ value, color: decodeColor(color) ?? FALLBACK_COLOR })) };
    }
    return { name: 'list', params: { kind: 'set', colors: [] } };
}

/** Colors for amino acid groups, based on Clustal (https://www.jalview.org/help/html/colourSchemes/clustal.html) */
const AminoGroupColors = {
    aromatic: decodeColor('#15A4A4')!,
    hydrophobic: decodeColor('#80A0F0')!,
    polar: decodeColor('#15C015')!,
    positive: decodeColor('#F01505')!,
    negative: decodeColor('#C048C0')!,
    proline: decodeColor('#C0C000')!,
    cysteine: decodeColor('#F08080')!,
    glycine: decodeColor('#F09048')!,
};

/** Colors for individual amino acids, based on Clustal (https://www.jalview.org/help/html/colourSchemes/clustal.html) */
const ResiduePropertyColors = {
    ...ResidueNameColors,
    HIS: AminoGroupColors.aromatic,
    TYR: AminoGroupColors.aromatic,
    ALA: AminoGroupColors.hydrophobic,
    VAL: AminoGroupColors.hydrophobic,
    LEU: AminoGroupColors.hydrophobic,
    ILE: AminoGroupColors.hydrophobic,
    MET: AminoGroupColors.hydrophobic,
    PHE: AminoGroupColors.hydrophobic,
    TRP: AminoGroupColors.hydrophobic,
    SER: AminoGroupColors.polar,
    THR: AminoGroupColors.polar,
    ASN: AminoGroupColors.polar,
    GLN: AminoGroupColors.polar,
    LYS: AminoGroupColors.positive,
    ARG: AminoGroupColors.positive,
    ASP: AminoGroupColors.negative,
    GLU: AminoGroupColors.negative,
    PRO: AminoGroupColors.proline,
    CYS: AminoGroupColors.cysteine,
    GLY: AminoGroupColors.glycine,
};

const MvsNamedColorListToMolstarName: Record<ColorListNameT, ColorListName> = {
    Blues: 'blues', Greens: 'greens', Greys: 'greys', Oranges: 'oranges', Purples: 'purples', Reds: 'reds',
    BuGn: 'blue-green', BuPu: 'blue-purple', GnBu: 'green-blue', OrRd: 'orange-red', PuBuGn: 'purple-blue-green',
    PuBu: 'purple-blue', PuRd: 'purple-red', RdPu: 'red-purple', YlGnBu: 'yellow-green-blue', YlGn: 'yellow-green',
    YlOrBr: 'yellow-orange-brown', YlOrRd: 'yellow-orange-red',
    BrBG: 'brown-white-green', PRGn: 'purple-green', PiYG: 'pink-yellow-green', PuOr: 'purple-orange', RdBu: 'red-blue',
    RdGy: 'red-grey', RdYlBu: 'red-yellow-blue', RdYlGn: 'red-yellow-green', Spectral: 'spectral',
    Cividis: 'cividis', Viridis: 'viridis', Inferno: 'inferno', Magma: 'magma', Plasma: 'plasma',
    Warm: 'warm', Cool: 'cool', CubehelixDefault: 'cubehelix-default', Turbo: 'turbo',
    Rainbow: 'rainbow', Sinebow: 'sinebow',
    Category10: 'category-10', Observable10: 'observable-10', Tableau10: 'tableau-10',
    Set1: 'set-1', Set2: 'set-2', Set3: 'set-3', Pastel1: 'pastel-1', Pastel2: 'pastel-2', Dark2: 'dark-2', Paired: 'paired', Accent: 'accent',
    Chainbow: 'turbo-no-black',
};

const MvsNamedColorDicts: Record<ColorDictNameT, Record<string, Color>> = {
    ElementSymbol: omitObjectKeys(ElementSymbolColors, ['C']),
    ResidueName: ResidueNameColors,
    ResidueProperties: ResiduePropertyColors,
};

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
    const params = node.params;
    switch (params.type) {
        case 'isosurface':
            const isoValue = typeof params.absolute_isovalue === 'number' ? Volume.IsoValue.absolute(params.absolute_isovalue) : Volume.IsoValue.relative(params.relative_isovalue ?? 0);
            const visuals: ('wireframe' | 'solid')[] = [];
            if (params.show_wireframe) visuals.push('wireframe');
            if (params.show_faces) visuals.push('solid');
            return {
                type: { name: 'isosurface', params: { alpha, isoValue, visuals } },
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