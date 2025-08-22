/**
 * Copyright (c) 2024-2025 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Adam Midlik <midlik@gmail.com>
 */

import { BaseGeometry } from '../../../mol-geo/geometry/base';
import { Lines } from '../../../mol-geo/geometry/lines/lines';
import { LinesBuilder } from '../../../mol-geo/geometry/lines/lines-builder';
import { addFixedCountDashedCylinder, addSimpleCylinder, BasicCylinderProps } from '../../../mol-geo/geometry/mesh/builder/cylinder';
import { addEllipsoid } from '../../../mol-geo/geometry/mesh/builder/ellipsoid';
import { Mesh } from '../../../mol-geo/geometry/mesh/mesh';
import { MeshBuilder } from '../../../mol-geo/geometry/mesh/mesh-builder';
import { Text } from '../../../mol-geo/geometry/text/text';
import { TextBuilder } from '../../../mol-geo/geometry/text/text-builder';
import { Box, BoxCage } from '../../../mol-geo/primitive/box';
import { Circle } from '../../../mol-geo/primitive/circle';
import { Primitive } from '../../../mol-geo/primitive/primitive';
import { StringLike } from '../../../mol-io/common/string-like';
import { Box3D, Sphere3D } from '../../../mol-math/geometry';
import { Mat4, Vec3 } from '../../../mol-math/linear-algebra';
import { radToDeg } from '../../../mol-math/misc';
import { Shape } from '../../../mol-model/shape';
import { Structure, StructureElement, StructureSelection } from '../../../mol-model/structure';
import { StructureQueryHelper } from '../../../mol-plugin-state/helpers/structure-query';
import { PluginStateObject as SO } from '../../../mol-plugin-state/objects';
import { PluginContext } from '../../../mol-plugin/context';
import { ShapeRepresentation } from '../../../mol-repr/shape/representation';
import { Expression } from '../../../mol-script/language/expression';
import { StateObject, StateTransformer } from '../../../mol-state';
import { Task } from '../../../mol-task';
import { round } from '../../../mol-util';
import { range } from '../../../mol-util/array';
import { Asset } from '../../../mol-util/assets';
import { Color } from '../../../mol-util/color';
import { MarkerActions } from '../../../mol-util/marker-action';
import { ParamDefinition as PD } from '../../../mol-util/param-definition';
import { capitalize } from '../../../mol-util/string';
import { rowsToExpression, rowToExpression } from '../helpers/selections';
import { collectMVSReferences, decodeColor, isDefined } from '../helpers/utils';
import { addParamDefaults } from '../tree/generic/params-schema';
import { MolstarNode, MolstarNodeParams, MolstarSubtree } from '../tree/molstar/molstar-tree';
import { MVSNode, MVSTreeSchema } from '../tree/mvs/mvs-tree';
import { isComponentExpression, isPrimitiveComponentExpressions, isVector3, PrimitivePositionT } from '../tree/mvs/param-types';
import { MVSTransform } from './annotation-structure-component';


type PrimitivesParams = MolstarNode<'primitives'>['params']

type _PrimitiveParams = MolstarNode<'primitive'>['params']
type PrimitiveKind = _PrimitiveParams['kind']
type PrimitiveParams<T extends PrimitiveKind = PrimitiveKind> = Extract<_PrimitiveParams, { kind: T }>

export function getPrimitiveStructureRefs(primitives: MolstarSubtree<'primitives'>) {
    const refs = new Set<string>();
    for (const c of primitives.children ?? []) {
        if (c.kind !== 'primitive') continue;
        const p = c.params;
        Builders[p.kind].resolveRefs?.(p, refs);
    }
    return refs;
}

export class MVSPrimitivesData extends SO.Create<PrimitiveBuilderContext>({ name: 'Primitive Data', typeClass: 'Object' }) { }
export class MVSPrimitiveShapes extends SO.Create<{ mesh?: Shape<Mesh>, labels?: Shape<Text> }>({ name: 'Primitive Shapes', typeClass: 'Object' }) { }

export type MVSDownloadPrimitiveData = typeof MVSDownloadPrimitiveData
export const MVSDownloadPrimitiveData = MVSTransform({
    name: 'mvs-download-primitive-data',
    display: { name: 'MVS Primitives' },
    from: [SO.Root, SO.Molecule.Structure],
    to: MVSPrimitivesData,
    params: {
        uri: PD.Url('', { isHidden: true }),
        format: PD.Text<'mvs-node-json'>('mvs-node-json', { isHidden: true })
    },
})({
    apply({ a, params, cache }, plugin: PluginContext) {
        return Task.create('Download Primitive Data', async ctx => {
            const url = Asset.getUrlAsset(plugin.managers.asset, params.uri);
            const asset = await plugin.managers.asset.resolve(url, 'string').runInContext(ctx);
            const node = JSON.parse(StringLike.toString(asset.data)) as MolstarSubtree<'primitives'>;
            (cache as any).asset = asset;
            return new MVSPrimitivesData({
                node,
                defaultStructure: SO.Molecule.Structure.is(a) ? a.data : undefined,
                structureRefs: {},
                primitives: getPrimitives(node),
                options: { ...node.params },
                positionCache: new Map(),
                instances: getInstances(node.params),
            }, { label: 'Primitive Data' });
        });
    },
    dispose({ cache }) {
        ((cache as any)?.asset as Asset.Wrapper | undefined)?.dispose();
    },
});

/* Cannot use MolstarSubtree<'primitives'>> because information about type of children would be lost and cause TypeScript errors in dependent code */
interface PrimitivesSubtree {
    kind: 'primitives',
    params: MolstarNodeParams<'primitives'>,
    children?: {
        kind: 'primitive',
        params: MolstarNodeParams<'primitive'>,
    }[],
}

export type MVSInlinePrimitiveData = typeof MVSInlinePrimitiveData
export const MVSInlinePrimitiveData = MVSTransform({
    name: 'mvs-inline-primitive-data',
    display: { name: 'MVS Primitives' },
    from: [SO.Root, SO.Molecule.Structure],
    to: MVSPrimitivesData,
    params: {
        node: PD.Value<PrimitivesSubtree>({
            kind: 'primitives',
            params: addParamDefaults(MVSTreeSchema.nodes.primitives.params, {}),
        }, { isHidden: true }),
    },
})({
    apply({ a, params }) {
        return new MVSPrimitivesData({
            node: params.node,
            defaultStructure: SO.Molecule.Structure.is(a) ? a.data : undefined,
            structureRefs: {},
            primitives: getPrimitives(params.node),
            options: { ...params.node.params },
            positionCache: new Map(),
            instances: getInstances(params.node.params),
        }, { label: 'Primitive Data' });
    }
});

export type MVSBuildPrimitiveShape = typeof MVSBuildPrimitiveShape
export const MVSBuildPrimitiveShape = MVSTransform({
    name: 'mvs-build-primitive-shape',
    display: { name: 'MVS Primitives' },
    from: MVSPrimitivesData,
    to: SO.Shape.Provider,
    params: {
        kind: PD.Text<'mesh' | 'labels' | 'lines'>('mesh')
    }
})({
    apply({ a, params, dependencies }) {
        const structureRefs = dependencies ? collectMVSReferences([SO.Molecule.Structure], dependencies) : {};
        const context: PrimitiveBuilderContext = { ...a.data, structureRefs };

        const snapshotKey = { snapshotKey: { ...SnapshotKey, defaultValue: a.data.options?.snapshot_key ?? '' } };
        const markdownCommands = { markdownCommands: { ...MarkdownCommands, defaultValue: a.data.node?.custom?.molstar_markdown_commands } };

        const label = capitalize(params.kind);
        if (params.kind === 'mesh') {
            if (!hasPrimitiveKind(a.data, 'mesh')) return StateObject.Null;

            const customMeshParams = a.data.node.custom?.molstar_mesh_params;
            return new SO.Shape.Provider({
                label,
                data: context,
                params: {
                    ...PD.withDefaults(Mesh.Params, { alpha: a.data.options?.opacity ?? 1, ...customMeshParams }),
                    ...snapshotKey,
                    ...markdownCommands,
                },
                getShape: (_, data, __, prev: any) => buildPrimitiveMesh(data, prev?.geometry),
                geometryUtils: Mesh.Utils,
            }, { label });
        } else if (params.kind === 'labels') {
            if (!hasPrimitiveKind(a.data, 'label')) return StateObject.Null;

            const options = a.data.options;
            const bgColor = options?.label_background_color;
            const customLabelParams = a.data.node.custom?.molstar_label_params;
            return new SO.Shape.Provider({
                label,
                data: context,
                params: {
                    ...PD.withDefaults(DefaultLabelParams, {
                        alpha: a.data.options?.label_opacity ?? 1,
                        attachment: options?.label_attachment ?? 'middle-center',
                        tether: options?.label_show_tether ?? false,
                        tetherLength: options?.label_tether_length ?? 1,
                        background: isDefined(bgColor),
                        backgroundColor: isDefined(bgColor) ? decodeColor(bgColor) : undefined,
                        ...customLabelParams,
                    }),
                    ...snapshotKey,
                    ...markdownCommands,
                },
                getShape: (_, data, props, prev: any) => buildPrimitiveLabels(data, prev?.geometry, props),
                geometryUtils: Text.Utils,
            }, { label });
        } else if (params.kind === 'lines') {
            if (!hasPrimitiveKind(a.data, 'line')) return StateObject.Null;

            const customLineParams = a.data.node.custom?.molstar_line_params;
            return new SO.Shape.Provider({
                label,
                data: context,
                params: {
                    ...PD.withDefaults(Lines.Params, { alpha: a.data.options?.opacity ?? 1, ...customLineParams }),
                    ...snapshotKey,
                    ...markdownCommands,
                },
                getShape: (_, data, __, prev: any) => buildPrimitiveLines(data, prev?.geometry),
                geometryUtils: Lines.Utils,
            }, { label });
        }

        return StateObject.Null;
    }
});

export const MVSShapeRepresentation3D = MVSTransform({
    name: 'shape-representation-3d',
    display: '3D Representation',
    from: SO.Shape.Provider,
    to: SO.Shape.Representation3D,
    params: (a, ctx: PluginContext) => {
        return a ? a.data.params : BaseGeometry.Params;
    }
})({
    canAutoUpdate() {
        return true;
    },
    apply({ a, params }) {
        return Task.create('Shape Representation', async ctx => {
            const props = { ...PD.getDefaultValues(a.data.params), ...params };
            const repr = ShapeRepresentation(a.data.getShape, a.data.geometryUtils);
            await repr.createOrUpdate(props, a.data.data).runInContext(ctx);

            const pickable = !!(params as any).snapshotKey?.trim() || !!(params as any).markdownCommands;
            if (pickable) {
                repr.setState({ pickable, markerActions: MarkerActions.Highlighting });
            }

            return new SO.Shape.Representation3D({ repr, sourceData: a.data }, { label: a.data.label });
        });
    },
    update({ a, b, newParams }) {
        return Task.create('Shape Representation', async ctx => {
            const props = { ...b.data.repr.props, ...newParams };
            await b.data.repr.createOrUpdate(props, a.data.data).runInContext(ctx);
            b.data.sourceData = a.data;

            const pickable = !!(newParams as any).snapshotKey?.trim() || !!(newParams as any).markdownCommands;
            if (pickable) {
                b.data.repr.setState({ pickable, markerActions: MarkerActions.Highlighting });
            }

            return StateTransformer.UpdateResult.Updated;
        });
    }
});

const SnapshotKey = PD.Text('', { isEssential: true, disableInteractiveUpdates: true, description: 'Activate the snapshot with the provided key when clicking on the label' });
const MarkdownCommands = PD.Value<any>(undefined, { isHidden: true });

/* **************************************************** */

class GroupManager {
    private current = -1;
    groupToNodeMap = new Map<number, MVSNode<'primitive'>>();
    sizes = new Map<number, number>();
    colors = new Map<number, number>();
    tooltips = new Map<number, string>();

    allocateSingle(node: MVSNode<'primitive'>) {
        const group = ++this.current;
        this.groupToNodeMap.set(group, node);
        return group;
    }

    allocateMany(node: MVSNode<'primitive'>, groups: number[]) {
        const newGroups = new Map<number, number>();
        const base = this.current;
        for (const g of groups) {
            if (newGroups.has(g)) continue;
            const group = base + newGroups.size + 1;
            this.groupToNodeMap.set(group, node);
            newGroups.set(g, group);
        }
        this.current += newGroups.size + 1;
        return newGroups;
    }

    updateColor(group: number, color?: string | null) {
        const c = decodeColor(color);
        if (typeof c === 'number') this.colors.set(group, c);
    }

    updateTooltip(group: number, tooltip?: string | null) {
        if (typeof tooltip === 'string') this.tooltips.set(group, tooltip);
    }

    updateSize(group: number, size?: number | null) {
        if (typeof size === 'number') this.sizes.set(group, size);
    }
}

interface PrimitiveBuilderContext {
    node: MolstarNode<'primitives'>;
    defaultStructure?: Structure;
    structureRefs: Record<string, Structure | undefined>;
    primitives: MolstarNode<'primitive'>[];
    options: PrimitivesParams;
    positionCache: Map<string, [isDefined: boolean, Sphere3D, Box3D]>;
    instances: Mat4[] | undefined;
    emptySelectionWarningPrinted?: boolean;
}

function printEmptySelectionWarning(ctx: PrimitiveBuilderContext, position: PrimitivePositionT): void {
    if (!ctx.emptySelectionWarningPrinted) {
        console.warn('Some primitives use positions which refer to empty substructure, not showing these primitives.', position, '(There may be more)');
        ctx.emptySelectionWarningPrinted = true;
    }
}

interface MeshBuilderState {
    groups: GroupManager;
    mesh: MeshBuilder.State;
}

interface LabelBuilderState {
    groups: GroupManager;
    labels: TextBuilder;
}

interface LineBuilderState {
    groups: GroupManager;
    lines: LinesBuilder;
}

const BaseLabelProps: PD.Values<Text.Params> = {
    ...PD.getDefaultValues(Text.Params),
    attachment: 'middle-center',
    fontQuality: 3,
    fontWeight: 'normal',
    borderWidth: 0.15,
    borderColor: Color(0x0),
    background: false,
    backgroundOpacity: 0.5,
    tether: false,
};
const DefaultLabelParams = PD.withDefaults(Text.Params, BaseLabelProps);

interface PrimitiveBuilder {
    builders: {
        mesh?: (context: PrimitiveBuilderContext, state: MeshBuilderState, node: MVSNode<'primitive'>, params: any) => void,
        line?: (context: PrimitiveBuilderContext, state: LineBuilderState, node: MVSNode<'primitive'>, params: any) => void,
        label?: (context: PrimitiveBuilderContext, state: LabelBuilderState, node: MVSNode<'primitive'>, params: any) => void,
    }
    isApplicable?: {
        mesh?: ((primitive: any, context: PrimitiveBuilderContext) => boolean),
        line?: ((primitive: any, context: PrimitiveBuilderContext) => boolean),
        label?: ((primitive: any, context: PrimitiveBuilderContext) => boolean),
    },
    resolveRefs?: (params: any, refs: Set<string>) => void,
}

const Builders: Record<PrimitiveParams['kind'], PrimitiveBuilder> = {
    mesh: {
        builders: {
            mesh: addMesh,
            line: addMeshWireframe,
        },
        isApplicable: {
            mesh: (m: PrimitiveParams<'mesh'>) => m.show_triangles,
            line: (m: PrimitiveParams<'mesh'>) => m.show_wireframe,
        },
    },
    lines: {
        builders: {
            line: addLines,
        },
    },
    tube: {
        builders: {
            mesh: addTubeMesh,
        },
        resolveRefs: resolveLineRefs,
    },
    arrow: {
        builders: {
            mesh: addArrowMesh,
        },
        resolveRefs: (params: PrimitiveParams<'arrow'>, refs: Set<string>) => {
            addRef(params.start, refs);
            if (params.end) addRef(params.end, refs);
        },
    },
    label: {
        builders: {
            label: addPrimitiveLabel,
        },
        resolveRefs: resolveLabelRefs,
    },
    distance_measurement: {
        builders: {
            mesh: addDistanceMesh,
            label: addDistanceLabel,
        },
        resolveRefs: resolveLineRefs,
    },
    angle_measurement: {
        builders: {
            mesh: addAngleMesh,
            label: addAngleLabel,
        },
        resolveRefs: (params: PrimitiveParams<'angle_measurement'>, refs: Set<string>) => {
            addRef(params.a, refs);
            addRef(params.b, refs);
            addRef(params.c, refs);
        },
    },
    ellipse: {
        builders: {
            mesh: addEllipseMesh,
        },
        resolveRefs: (params: PrimitiveParams<'ellipse'>, refs: Set<string>) => {
            addRef(params.center, refs);
            if (params.major_axis_endpoint) addRef(params.major_axis_endpoint, refs);
            if (params.minor_axis_endpoint) addRef(params.minor_axis_endpoint, refs);
        },
    },
    ellipsoid: {
        builders: {
            mesh: addEllipsoidMesh,
        },
        resolveRefs: (params: PrimitiveParams<'ellipsoid'>, refs: Set<string>) => {
            addRef(params.center, refs);
            if (params.major_axis_endpoint) addRef(params.major_axis_endpoint, refs);
            if (params.minor_axis_endpoint) addRef(params.minor_axis_endpoint, refs);
        },
    },
    box: {
        builders: {
            mesh: addBoxMesh,
        },
        resolveRefs: (params: PrimitiveParams<'box'>, refs: Set<string>) => {
            addRef(params.center, refs);
        },
    }
};


function getPrimitives(primitives: MolstarSubtree<'primitives'>) {
    return (primitives.children ?? []).filter(c => c.kind === 'primitive') as MolstarNode<'primitive'>[];
}

function addRef(position: PrimitivePositionT, refs: Set<string>) {
    if (isPrimitiveComponentExpressions(position) && position.structure_ref) {
        refs.add(position.structure_ref);
    }
}

function hasPrimitiveKind(context: PrimitiveBuilderContext, kind: 'mesh' | 'line' | 'label') {
    for (const c of context.primitives) {
        const params = c.params;
        const b = Builders[params.kind];
        const builderFunction = b.builders[kind];
        if (builderFunction) {
            const test = b.isApplicable?.[kind];
            if (test === undefined || test(params, context)) {
                return true;
            }
        }
    }
    return false;
}

/** Save resolved position into `targetPosition`.
 * Return `true` if the resolved position is defined (i.e. vector or non-empty selection);
 * return `false` if the resolved position is not defined (i.e. empty selection). */
function resolveBasePosition(context: PrimitiveBuilderContext, position: PrimitivePositionT, targetPosition: Vec3): boolean {
    return resolvePosition(context, position, targetPosition, undefined, undefined);
}

const _EmptySphere = Sphere3D.zero();
const _EmptyBox = Box3D.zero();

/** Save resolved position into `targetPosition`, `targetSphere`, `targetBox`.
 * Return `true` if the resolved position is defined (i.e. vector or non-empty selection);
 * return `false` if the resolved position is not defined (i.e. empty selection). */
function resolvePosition(context: PrimitiveBuilderContext, position: PrimitivePositionT, targetPosition: Vec3 | undefined, targetSphere: Sphere3D | undefined, targetBox: Box3D | undefined): boolean {
    let expr: Expression | undefined;
    let pivotRef: string | undefined;

    if (isVector3(position)) {
        if (targetPosition) Vec3.copy(targetPosition, position as any);
        if (targetSphere) Sphere3D.set(targetSphere, position as any, 0);
        if (targetBox) Box3D.set(targetBox, position as any, position as any);
        return true;
    }

    if (isPrimitiveComponentExpressions(position)) {
        // TODO: take schema into account for possible optimization
        expr = rowsToExpression(position.expressions!);
        pivotRef = position.structure_ref;
    } else if (isComponentExpression(position)) {
        expr = rowToExpression(position);
    }

    if (!expr) {
        console.error('Invalid expression', position);
        throw new Error('Invalid primitive potition expression, see console for details.');
    }

    const pivot = !pivotRef ? context.defaultStructure : context.structureRefs[pivotRef];
    if (!pivot) {
        throw new Error(`Structure with ref '${pivotRef ?? '<default>'}' not found.`);
    }

    const cacheKey = JSON.stringify(position);
    if (context.positionCache.has(cacheKey)) {
        const [isDefined, sphere, box] = context.positionCache.get(cacheKey)!;
        if (targetPosition) Vec3.copy(targetPosition, sphere.center);
        if (targetSphere) Sphere3D.copy(targetSphere, sphere);
        if (targetBox) Box3D.copy(targetBox, box);
        return isDefined;
    }

    const { selection } = StructureQueryHelper.createAndRun(pivot, expr);

    let box: Box3D;
    let sphere: Sphere3D;
    let isDefined: boolean;

    if (StructureSelection.isEmpty(selection)) {
        if (targetPosition) Vec3.set(targetPosition, 0, 0, 0);
        box = _EmptyBox;
        sphere = _EmptySphere;
        isDefined = false;
        printEmptySelectionWarning(context, position);
    } else {
        const loci = StructureSelection.toLociWithSourceUnits(selection);
        const boundary = StructureElement.Loci.getBoundary(loci);
        if (targetPosition) Vec3.copy(targetPosition, boundary.sphere.center);
        box = boundary.box;
        sphere = boundary.sphere;
        isDefined = true;
    }

    if (targetSphere) Sphere3D.copy(targetSphere, sphere);
    if (targetBox) Box3D.copy(targetBox, box);

    context.positionCache.set(cacheKey, [isDefined, sphere, box]);
    return isDefined;
}

function getInstances(options: PrimitivesParams | undefined): Mat4[] | undefined {
    if (!options?.instances?.length) return undefined;
    return options.instances.map(i => Mat4.fromArray(Mat4(), i, 0));
}

function buildPrimitiveMesh(context: PrimitiveBuilderContext, prev?: Mesh): Shape<Mesh> {
    const meshBuilder = MeshBuilder.createState(1024, 1024, prev);
    const state: MeshBuilderState = { groups: new GroupManager(), mesh: meshBuilder };

    meshBuilder.currentGroup = -1;

    for (const c of context.primitives) {
        const p = c.params;
        const b = Builders[p.kind];
        if (!b) {
            console.warn(`Primitive ${p.kind} not supported`);
            continue;
        }
        b.builders.mesh?.(context, state, c, p);
    }

    const { colors, tooltips } = state.groups;
    const tooltip = context.options?.tooltip ?? '';
    const color = decodeColor(context.options?.color) ?? Color(0);

    return Shape.create(
        'Mesh',
        {
            kind: 'mvs-primitives',
            node: context.node,
            groupToNode: state.groups.groupToNodeMap,
        },
        MeshBuilder.getMesh(meshBuilder),
        (g) => colors.get(g) as Color ?? color,
        (g) => 1,
        (g) => tooltips.get(g) ?? tooltip,
        context.instances,
    );
}

function buildPrimitiveLines(context: PrimitiveBuilderContext, prev?: Lines): Shape<Lines> {
    const linesBuilder = LinesBuilder.create(1024, 1024, prev);
    const state: LineBuilderState = { groups: new GroupManager(), lines: linesBuilder };

    for (const c of context.primitives) {
        const p = c.params;
        const b = Builders[p.kind];
        if (!b) {
            console.warn(`Primitive ${p.kind} not supported`);
            continue;
        }
        b.builders.line?.(context, state, c, p);
    }

    const { colors, sizes, tooltips } = state.groups;
    const tooltip = context.options?.tooltip ?? '';
    const color = decodeColor(context.options?.color) ?? Color(0);

    return Shape.create(
        'Lines',
        {
            kind: 'mvs-primitives',
            node: context.node,
            groupToNode: state.groups.groupToNodeMap,
        },
        linesBuilder.getLines(),
        (g) => colors.get(g) as Color ?? color,
        (g) => sizes.get(g) ?? 1,
        (g) => tooltips.get(g) ?? tooltip,
        context.instances,
    );
}

function buildPrimitiveLabels(context: PrimitiveBuilderContext, prev: Text | undefined, props: PD.Values<Text.Params>): Shape<Text> {
    const labelsBuilder = TextBuilder.create({
        ...BaseLabelProps,
        ...props,
    }, 1024, 1024, prev);
    const state: LabelBuilderState = { groups: new GroupManager(), labels: labelsBuilder };

    for (const c of context.primitives) {
        const p = c.params;
        const b = Builders[p.kind];
        if (!b) {
            console.warn(`Primitive ${p.kind} not supported`);
            continue;
        }
        b.builders.label?.(context, state, c, p);
    }

    const color = decodeColor(context.options?.label_color) ?? Color(0);
    const { colors, sizes, tooltips } = state.groups;

    return Shape.create(
        'Labels',
        {
            kind: 'mvs-primitives',
            node: context.node,
            groupToNode: state.groups.groupToNodeMap,
        },
        labelsBuilder.getText(),
        (g) => colors.get(g) as Color ?? color,
        (g) => sizes.get(g) ?? 1,
        (g) => tooltips.get(g) ?? '',
        context.instances,
    );
}

function addMeshFaces(context: PrimitiveBuilderContext, groups: GroupManager, node: MVSNode<'primitive'>, params: PrimitiveParams<'mesh'>, addFace: (mvsGroup: number, builderGroup: number, a: Vec3, b: Vec3, c: Vec3) => void) {
    const a = Vec3.zero();
    const b = Vec3.zero();
    const c = Vec3.zero();

    let { indices, vertices, triangle_groups } = params;
    const nTriangles = Math.floor(indices.length / 3);
    triangle_groups ??= range(nTriangles); // implicit grouping (triangle i = group i)
    const groupSet = groups.allocateMany(node, triangle_groups);

    for (let i = 0; i < nTriangles; i++) {
        const mvsGroup = triangle_groups[i];
        const builderGroup = groupSet.get(mvsGroup)!;
        Vec3.fromArray(a, vertices, 3 * indices[3 * i]);
        Vec3.fromArray(b, vertices, 3 * indices[3 * i + 1]);
        Vec3.fromArray(c, vertices, 3 * indices[3 * i + 2]);

        addFace(mvsGroup, builderGroup, a, b, c);
    }
}

function addMesh(context: PrimitiveBuilderContext, { groups, mesh }: MeshBuilderState, node: MVSNode<'primitive'>, params: PrimitiveParams<'mesh'>) {
    if (!params.show_triangles) return;

    const { group_colors, group_tooltips, color, tooltip } = params;

    addMeshFaces(context, groups, node, params, (mvsGroup, builderGroup, a, b, c) => {
        groups.updateColor(builderGroup, group_colors[mvsGroup] ?? color);
        groups.updateTooltip(builderGroup, group_tooltips[mvsGroup] ?? tooltip);
        mesh.currentGroup = builderGroup;
        MeshBuilder.addTriangle(mesh, a, b, c);
    });
    // this could be slightly improved by only updating color and tooltip once per group instead of once per triangle
}

function addMeshWireframe(context: PrimitiveBuilderContext, { groups, lines }: LineBuilderState, node: MVSNode<'primitive'>, params: PrimitiveParams<'mesh'>) {
    if (!params.show_wireframe) return;
    const width = params.wireframe_width;

    const { group_colors, group_tooltips, wireframe_color, color, tooltip } = params;

    addMeshFaces(context, groups, node, params, (mvsGroup, builderGroup, a, b, c) => {
        groups.updateColor(builderGroup, wireframe_color ?? group_colors[mvsGroup] ?? color);
        groups.updateTooltip(builderGroup, group_tooltips[mvsGroup] ?? tooltip);
        groups.updateSize(builderGroup, width);
        lines.add(a[0], a[1], a[2], b[0], b[1], b[2], builderGroup);
        lines.add(b[0], b[1], b[2], c[0], c[1], c[2], builderGroup);
        lines.add(c[0], c[1], c[2], a[0], a[1], a[2], builderGroup);
    });
}

function addLines(context: PrimitiveBuilderContext, { groups, lines }: LineBuilderState, node: MVSNode<'primitive'>, params: PrimitiveParams<'lines'>) {
    const a = Vec3.zero();
    const b = Vec3.zero();

    let { indices, vertices, line_groups, group_colors, group_tooltips, group_widths } = params;
    const width = params.width;

    const nLines = Math.floor(indices.length / 2);
    line_groups ??= range(nLines); // implicit grouping (line i = group i)
    const groupSet = groups.allocateMany(node, line_groups);

    for (let i = 0; i < nLines; i++) {
        const mvsGroup = line_groups[i];
        const builderGroup = groupSet.get(mvsGroup)!;
        groups.updateColor(builderGroup, group_colors[mvsGroup] ?? params.color);
        groups.updateTooltip(builderGroup, group_tooltips[mvsGroup] ?? params.tooltip);
        groups.updateSize(builderGroup, group_widths[mvsGroup] ?? width);

        Vec3.fromArray(a, vertices, 3 * indices[2 * i]);
        Vec3.fromArray(b, vertices, 3 * indices[2 * i + 1]);
        lines.add(a[0], a[1], a[2], b[0], b[1], b[2], builderGroup);
    }
}

function resolveLineRefs(params: PrimitiveParams<'tube' | 'distance_measurement'>, refs: Set<string>) {
    addRef(params.start, refs);
    addRef(params.end, refs);
}

const lStart = Vec3.zero();
const lEnd = Vec3.zero();

function addTubeMesh(context: PrimitiveBuilderContext, { groups, mesh }: MeshBuilderState, node: MVSNode<'primitive'>, params: PrimitiveParams<'tube'>, options?: { skipResolvePosition?: boolean }) {
    if (!options?.skipResolvePosition) {
        const startDefined = resolveBasePosition(context, params.start, lStart);
        const endDefined = resolveBasePosition(context, params.end, lEnd);
        if (!startDefined || !endDefined) return;
    }
    const radius = params.radius;

    const cylinderProps: BasicCylinderProps = {
        radiusBottom: radius,
        radiusTop: radius,
        topCap: true,
        bottomCap: true,
    };

    mesh.currentGroup = groups.allocateSingle(node);
    groups.updateColor(mesh.currentGroup, params.color);
    groups.updateTooltip(mesh.currentGroup, params.tooltip);

    if (params.dash_length) {
        const dist = Vec3.distance(lStart, lEnd);
        const count = Math.ceil(dist / (2 * params.dash_length));
        addFixedCountDashedCylinder(mesh, lStart, lEnd, 1.0, count, true, cylinderProps);
    } else {
        addSimpleCylinder(mesh, lStart, lEnd, cylinderProps);
    }
}

const ArrowState = {
    start: Vec3.zero(),
    end: Vec3.zero(),
    dir: Vec3.zero(),
    startCap: Vec3.zero(),
    endCap: Vec3.zero(),
};

function addArrowMesh(context: PrimitiveBuilderContext, { groups, mesh }: MeshBuilderState, node: MVSNode<'primitive'>, params: PrimitiveParams<'arrow'>) {
    const startDefined = resolveBasePosition(context, params.start, ArrowState.start);
    if (!startDefined) return;

    if (params.end) {
        const endDefined = resolveBasePosition(context, params.end, ArrowState.end);
        if (!endDefined) return;
    } else if (params.direction) {
        Vec3.add(ArrowState.end, ArrowState.start, params.direction as any as Vec3);
    } else {
        console.warn(`Primitive arrow does not contain "end" nor "distance". Not showing.`);
        return;
    }

    Vec3.sub(ArrowState.dir, ArrowState.end, ArrowState.start);
    Vec3.normalize(ArrowState.dir, ArrowState.dir);

    if (params.length) {
        Vec3.scaleAndAdd(ArrowState.end, ArrowState.start, ArrowState.dir, params.length);
    }

    const length = Vec3.distance(ArrowState.start, ArrowState.end);
    if (length < 1e-3) return;

    const tubeRadius = params.tube_radius;
    const tubeProps: BasicCylinderProps = {
        radiusBottom: tubeRadius,
        radiusTop: tubeRadius,
        topCap: !params.show_end_cap,
        bottomCap: !params.show_start_cap,
    };

    mesh.currentGroup = groups.allocateSingle(node);
    groups.updateColor(mesh.currentGroup, params.color);
    groups.updateTooltip(mesh.currentGroup, params.tooltip);

    if (params.show_start_cap) {
        const startRadius = params.start_cap_radius ?? 2 * tubeRadius;
        const startCapLength = params.start_cap_length ?? 2 * startRadius;
        Vec3.scaleAndAdd(ArrowState.startCap, ArrowState.start, ArrowState.dir, startCapLength);
        addSimpleCylinder(mesh, ArrowState.startCap, ArrowState.start, {
            radiusBottom: startRadius,
            radiusTop: 0,
            topCap: false,
            bottomCap: true,
            radialSegments: 12,
        });
    } else {
        Vec3.copy(ArrowState.startCap, ArrowState.start);
    }

    if (params.show_end_cap) {
        const endRadius = params.end_cap_radius ?? 2 * tubeRadius;
        const endCapLength = params.end_cap_length ?? 2 * endRadius;
        Vec3.scaleAndAdd(ArrowState.endCap, ArrowState.end, ArrowState.dir, -endCapLength);
        addSimpleCylinder(mesh, ArrowState.endCap, ArrowState.end, {
            radiusBottom: endRadius,
            radiusTop: 0,
            topCap: false,
            bottomCap: true,
            radialSegments: 12,
        });
    } else {
        Vec3.copy(ArrowState.endCap, ArrowState.end);
    }

    if (params.show_tube) {
        if (params.tube_dash_length) {
            const dist = Vec3.distance(ArrowState.startCap, ArrowState.endCap);
            const count = Math.ceil(dist / (2 * params.tube_dash_length));
            addFixedCountDashedCylinder(mesh, ArrowState.startCap, ArrowState.endCap, 1.0, count, true, tubeProps);
        } else {
            addSimpleCylinder(mesh, ArrowState.startCap, ArrowState.endCap, tubeProps);
        }
    }
}


/** Return distance in angstroms, or `undefined` if any of the endpoints corresponds to empty substructure.
 * This function also sets `lStart`, `lEnd` globals. */
function computeDistance(context: PrimitiveBuilderContext, start: PrimitivePositionT, end: PrimitivePositionT): number | undefined {
    const startDefined = resolveBasePosition(context, start, lStart);
    const endDefined = resolveBasePosition(context, end, lEnd);
    if (startDefined && endDefined) return Vec3.distance(lStart, lEnd);
    else return undefined;
}

// /** Return text for distance measurement label/tooltip. */
function distanceLabel(distance: number, params: PrimitiveParams<'distance_measurement'>): string {
    const distStr = `${round(distance, 2)} Å`;
    if (typeof params.label_template === 'string') return params.label_template.replace('{{distance}}', distStr);
    else return distStr;
}

function addDistanceMesh(context: PrimitiveBuilderContext, state: MeshBuilderState, node: MVSNode<'primitive'>, params: PrimitiveParams<'distance_measurement'>) {
    const distance = computeDistance(context, params.start, params.end); // sets lStart, lEnd
    if (distance === undefined) return; // empty substructure in measurement
    const tooltip = distanceLabel(distance, params);
    addTubeMesh(context, state, node, { ...params, tooltip } as any, { skipResolvePosition: true });
}

const labelPos = Vec3.zero();

function addDistanceLabel(context: PrimitiveBuilderContext, state: LabelBuilderState, node: MVSNode<'primitive'>, params: PrimitiveParams<'distance_measurement'>) {
    const { labels, groups } = state;
    const dist = computeDistance(context, params.start, params.end); // sets lStart, lEnd
    if (dist === undefined) return; // empty substructure in measurement

    let size: number | undefined;
    if (typeof params.label_size === 'number') {
        size = params.label_size;
    } else {
        size = Math.max(dist * (params.label_auto_size_scale), params.label_auto_size_min);
    }

    Vec3.add(labelPos, lStart, lEnd);
    Vec3.scale(labelPos, labelPos, 0.5);

    const group = groups.allocateSingle(node);
    groups.updateColor(group, params.label_color);
    groups.updateSize(group, size);

    labels.add(distanceLabel(dist, params), labelPos[0], labelPos[1], labelPos[2], 1.05 * (params.radius), 1, group);
}


const AngleState = {
    isDefined: false,
    a: Vec3(),
    b: Vec3(),
    c: Vec3(),
    ba: Vec3(),
    bc: Vec3(),
    labelPos: Vec3(),
    /** Sector radius */
    radius: 0,
    label: '',
};

function syncAngleState(context: PrimitiveBuilderContext, params: PrimitiveParams<'angle_measurement'>): void {
    const aDefined = resolveBasePosition(context, params.a, AngleState.a);
    const bDefined = resolveBasePosition(context, params.b, AngleState.b);
    const cDefined = resolveBasePosition(context, params.c, AngleState.c);
    AngleState.isDefined = aDefined && bDefined && cDefined;
    if (!AngleState.isDefined) return;

    Vec3.sub(AngleState.ba, AngleState.a, AngleState.b);
    Vec3.sub(AngleState.bc, AngleState.c, AngleState.b);
    const value = radToDeg(Vec3.angle(AngleState.ba, AngleState.bc));

    const angle = `${round(value, 2)}\u00B0`;
    AngleState.label = typeof params.label_template === 'string' ? params.label_template.replace('{{angle}}', angle) : angle;

    if (typeof params.section_radius === 'number') {
        AngleState.radius = params.section_radius;
    } else {
        AngleState.radius = Math.min(Vec3.magnitude(AngleState.ba), Vec3.magnitude(AngleState.bc));
        if (typeof params.section_radius_scale === 'number') {
            AngleState.radius *= params.section_radius_scale;
        }
    }
}

function addAngleMesh(context: PrimitiveBuilderContext, state: MeshBuilderState, node: MVSNode<'primitive'>, params: PrimitiveParams<'angle_measurement'>) {
    syncAngleState(context, params);
    if (!AngleState.isDefined) return; // empty substructure in measurement

    const { groups, mesh } = state;

    if (params.show_vector) {
        const radius = params.vector_radius ?? 0.05;
        const cylinderProps: BasicCylinderProps = {
            radiusBottom: radius,
            radiusTop: radius,
            topCap: true,
            bottomCap: true,
        };

        mesh.currentGroup = groups.allocateSingle(node);
        groups.updateColor(mesh.currentGroup, params.vector_color);
        groups.updateTooltip(mesh.currentGroup, AngleState.label);

        let count = Math.ceil(Vec3.magnitude(AngleState.ba) / (2 * radius));
        addFixedCountDashedCylinder(mesh, AngleState.a, AngleState.b, 1.0, count, true, cylinderProps);
        count = Math.ceil(Vec3.magnitude(AngleState.bc) / (2 * radius));
        addFixedCountDashedCylinder(mesh, AngleState.b, AngleState.c, 1.0, count, true, cylinderProps);
    }

    if (params.show_section) {
        const angle = Vec3.angle(AngleState.ba, AngleState.bc);
        Vec3.normalize(AngleState.ba, AngleState.ba);
        Vec3.normalize(AngleState.bc, AngleState.bc);
        Vec3.scale(AngleState.ba, AngleState.ba, AngleState.radius);
        Vec3.scale(AngleState.bc, AngleState.bc, AngleState.radius);

        addEllipseMesh(context, state, node, {
            kind: 'ellipse',
            as_circle: true,
            center: AngleState.b as any,
            major_axis_endpoint: null,
            major_axis: AngleState.ba as any,
            minor_axis_endpoint: null,
            minor_axis: AngleState.bc as any,
            radius_major: AngleState.radius,
            radius_minor: AngleState.radius,
            theta_start: 0,
            theta_end: angle,
            color: params.section_color,
            tooltip: AngleState.label,
        });
    }
}

function addAngleLabel(context: PrimitiveBuilderContext, state: LabelBuilderState, node: MVSNode<'primitive'>, params: PrimitiveParams<'angle_measurement'>) {
    const { labels, groups } = state;
    syncAngleState(context, params);
    if (!AngleState.isDefined) return; // empty substructure in measurement

    Vec3.normalize(AngleState.ba, AngleState.ba);
    Vec3.normalize(AngleState.bc, AngleState.bc);
    Vec3.scale(AngleState.ba, AngleState.ba, AngleState.radius);
    Vec3.scale(AngleState.bc, AngleState.bc, AngleState.radius);

    let size: number | undefined;
    if (typeof params.label_size === 'number') {
        size = params.label_size;
    } else {
        size = Math.max(AngleState.radius * (params.label_auto_size_scale), params.label_auto_size_min);
    }

    Vec3.add(AngleState.labelPos, AngleState.ba, AngleState.bc);
    Vec3.normalize(AngleState.labelPos, AngleState.labelPos);
    Vec3.scale(AngleState.labelPos, AngleState.labelPos, AngleState.radius);
    Vec3.add(AngleState.labelPos, AngleState.labelPos, AngleState.b);

    const group = groups.allocateSingle(node);
    groups.updateColor(group, params.label_color);
    groups.updateSize(group, size);

    labels.add(AngleState.label, AngleState.labelPos[0], AngleState.labelPos[1], AngleState.labelPos[2], 1, 1, group);
}

function resolveLabelRefs(params: PrimitiveParams<'label'>, refs: Set<string>) {
    addRef(params.position, refs);
}

const PrimitiveLabelState = {
    position: Vec3.zero(),
    sphere: Sphere3D.zero(),
};

function addPrimitiveLabel(context: PrimitiveBuilderContext, state: LabelBuilderState, node: MVSNode<'primitive'>, params: PrimitiveParams<'label'>) {
    const { labels, groups } = state;
    const positionDefined = resolvePosition(context, params.position, PrimitiveLabelState.position, PrimitiveLabelState.sphere, undefined);
    if (!positionDefined) return;

    const group = groups.allocateSingle(node);
    groups.updateColor(group, params.label_color);
    groups.updateSize(group, params.label_size);

    const offset = PrimitiveLabelState.sphere.radius + params.label_offset;
    labels.add(params.text, PrimitiveLabelState.position[0], PrimitiveLabelState.position[1], PrimitiveLabelState.position[2], offset, 1, group);
}

const circleCache = new Map<string, Primitive>();

function getCircle(options: { thetaStart?: number, thetaEnd?: number }) {
    const key = JSON.stringify(options);
    if (circleCache.has(key)) return circleCache.get(key)!;
    const thetaLength = (options.thetaEnd ?? 2 * Math.PI) - (options.thetaStart ?? 0);
    if (Math.abs(thetaLength) < 1e-3) return null;

    const circle = Circle({
        radius: 1,
        thetaStart: options.thetaStart ?? 0,
        thetaLength,
        segments: Math.ceil(2 * Math.PI / thetaLength * 64),
    });
    circleCache.set(key, circle);
    return circle;
}

const EllipseState = {
    centerPos: Vec3.zero(),
    majorPos: Vec3.zero(),
    minorPos: Vec3.zero(),
    majorAxis: Vec3.zero(),
    minorAxis: Vec3.zero(),
    scale: Vec3.zero(),
    normal: Vec3.zero(),
    scaleXform: Mat4.identity(),
    rotationXform: Mat4.identity(),
    translationXform: Mat4.identity(),
    xform: Mat4.identity(),
};


function addEllipseMesh(context: PrimitiveBuilderContext, state: MeshBuilderState, node: MVSNode<'primitive'>, params: PrimitiveParams<'ellipse'>) {
    // Unit circle in the XZ plane (Y up)
    // X = minor axis, Y = normal, Z = major axis

    const circle = getCircle({ thetaStart: params.theta_start, thetaEnd: params.theta_end });
    if (!circle) return;

    const centerDefined = resolvePosition(context, params.center, EllipseState.centerPos, undefined, undefined);
    if (!centerDefined) return;

    if (params.major_axis_endpoint) {
        const endpointDefined = resolvePosition(context, params.major_axis_endpoint, EllipseState.majorPos, undefined, undefined);
        if (!endpointDefined) return;
        Vec3.sub(EllipseState.majorAxis, EllipseState.majorPos, EllipseState.centerPos);
    } else {
        Vec3.copy(EllipseState.majorAxis, params.major_axis as any as Vec3);
    }

    if (params.minor_axis_endpoint) {
        const endpointDefined = resolvePosition(context, params.minor_axis_endpoint, EllipseState.minorPos, undefined, undefined);
        if (!endpointDefined) return;
        Vec3.sub(EllipseState.minorAxis, EllipseState.minorPos, EllipseState.centerPos);
    } else {
        Vec3.copy(EllipseState.minorAxis, params.minor_axis as any as Vec3);
    }

    const { mesh, groups } = state;

    // Translation
    Mat4.fromTranslation(EllipseState.translationXform, EllipseState.centerPos);

    // Scale
    if (params.as_circle) {
        const r = params.radius_major ?? Vec3.magnitude(EllipseState.majorAxis);
        Vec3.set(EllipseState.scale, r, 1, r);
    } else {
        const major = params.radius_major ?? Vec3.magnitude(EllipseState.majorAxis);
        const minor = params.radius_minor ?? Vec3.magnitude(EllipseState.minorAxis);
        Vec3.set(EllipseState.scale, minor, 1, major);
    }
    Mat4.fromScaling(EllipseState.scaleXform, EllipseState.scale);

    // Rotation
    Vec3.normalize(EllipseState.minorAxis, EllipseState.minorAxis);
    Vec3.normalize(EllipseState.majorAxis, EllipseState.majorAxis);
    Vec3.cross(EllipseState.normal, EllipseState.majorAxis, EllipseState.minorAxis);

    Mat4.targetTo(EllipseState.rotationXform, Vec3.origin, EllipseState.majorAxis, EllipseState.normal);
    Mat4.mul(EllipseState.rotationXform, EllipseState.rotationXform, Mat4.rotY180);

    // Final xform
    Mat4.mul3(EllipseState.xform, EllipseState.translationXform, EllipseState.rotationXform, EllipseState.scaleXform);

    mesh.currentGroup = groups.allocateSingle(node);
    groups.updateColor(mesh.currentGroup, params.color);
    groups.updateTooltip(mesh.currentGroup, params.tooltip);

    MeshBuilder.addPrimitive(mesh, EllipseState.xform, circle);
    MeshBuilder.addPrimitiveFlipped(mesh, EllipseState.xform, circle);
}

const EllipsoidState = {
    centerPos: Vec3.zero(),
    majorPos: Vec3.zero(),
    minorPos: Vec3.zero(),
    majorAxis: Vec3.zero(),
    minorAxis: Vec3.zero(),
    sphere: Sphere3D.zero(),
    radius: Vec3.zero(),
    extent: Vec3.zero(),
    up: Vec3.zero(),
};


function addEllipsoidMesh(context: PrimitiveBuilderContext, state: MeshBuilderState, node: MVSNode<'primitive'>, params: PrimitiveParams<'ellipsoid'>) {
    const centerDefined = resolvePosition(context, params.center, EllipsoidState.centerPos, EllipsoidState.sphere, undefined);
    if (!centerDefined) return;

    if (params.major_axis_endpoint) {
        const endpointDefined = resolvePosition(context, params.major_axis_endpoint, EllipsoidState.majorPos, undefined, undefined);
        if (!endpointDefined) return;
        Vec3.sub(EllipsoidState.majorAxis, EllipsoidState.majorPos, EllipsoidState.centerPos);
    } else if (params.major_axis) {
        Vec3.copy(EllipsoidState.majorAxis, params.major_axis as any as Vec3);
    } else {
        Vec3.copy(EllipsoidState.majorAxis, Vec3.unitX);
    }

    if (params.minor_axis_endpoint) {
        const endpointDefined = resolvePosition(context, params.minor_axis_endpoint, EllipsoidState.minorPos, undefined, undefined);
        if (!endpointDefined) return;
        Vec3.sub(EllipsoidState.minorAxis, EllipsoidState.minorPos, EllipsoidState.centerPos);
    } else if (params.minor_axis) {
        Vec3.copy(EllipsoidState.minorAxis, params.minor_axis as any as Vec3);
    } else {
        Vec3.copy(EllipsoidState.minorAxis, Vec3.unitY);
    }

    if (typeof params.radius === 'number') {
        Vec3.set(EllipsoidState.radius, params.radius, params.radius, params.radius);
    } else if (params.radius) {
        Vec3.copy(EllipsoidState.radius, params.radius as any as Vec3);
    } else {
        const r = EllipsoidState.sphere.radius;
        Vec3.set(EllipsoidState.radius, r, r, r);
    }

    if (typeof params.radius_extent === 'number') {
        Vec3.set(EllipsoidState.extent, params.radius_extent, params.radius_extent, params.radius_extent);
    } else if (params.radius_extent) {
        Vec3.copy(EllipsoidState.extent, params.radius_extent as any as Vec3);
    } else {
        Vec3.set(EllipsoidState.extent, 0, 0, 0);
    }

    Vec3.add(EllipsoidState.radius, EllipsoidState.radius, EllipsoidState.extent);

    const { mesh, groups } = state;

    mesh.currentGroup = groups.allocateSingle(node);
    groups.updateColor(mesh.currentGroup, params.color);
    groups.updateTooltip(mesh.currentGroup, params.tooltip);

    Vec3.normalize(EllipsoidState.majorAxis, EllipsoidState.majorAxis);
    Vec3.normalize(EllipsoidState.minorAxis, EllipsoidState.minorAxis);
    Vec3.cross(EllipsoidState.up, EllipsoidState.majorAxis, EllipsoidState.minorAxis);

    addEllipsoid(mesh, EllipsoidState.centerPos, EllipsoidState.up, EllipsoidState.minorAxis, EllipsoidState.radius, 3);
}


const BoxState = {
    center: Vec3.zero(),
    boundary: Box3D.zero(),
    size: Vec3.zero(),
    cage: BoxCage(),
    translationXform: Mat4.identity(),
    scaleXform: Mat4.identity(),
    xform: Mat4.identity(),
};

function addBoxMesh(context: PrimitiveBuilderContext, state: MeshBuilderState, node: MVSNode<'primitive'>, params: PrimitiveParams<'box'>) {
    if (!params.show_edges && !params.show_faces) return;

    const positionDefined = resolvePosition(context, params.center, BoxState.center, undefined, BoxState.boundary);
    if (!positionDefined) return;
    if (params.extent) {
        Box3D.expand(BoxState.boundary, BoxState.boundary, params.extent as unknown as Vec3);
    }

    if (Box3D.volume(BoxState.boundary) < 1e-3) return;

    const { mesh, groups } = state;

    Mat4.fromScaling(BoxState.scaleXform, Box3D.size(BoxState.size, BoxState.boundary));
    Mat4.fromTranslation(BoxState.translationXform, BoxState.center);
    Mat4.mul(BoxState.xform, BoxState.translationXform, BoxState.scaleXform);

    if (params.show_faces) {
        mesh.currentGroup = groups.allocateSingle(node);
        groups.updateColor(mesh.currentGroup, params.face_color);
        groups.updateTooltip(mesh.currentGroup, params.tooltip);
        MeshBuilder.addPrimitive(mesh, BoxState.xform, Box());
    }

    if (params.show_edges) {
        mesh.currentGroup = groups.allocateSingle(node);
        groups.updateColor(mesh.currentGroup, params.edge_color);
        groups.updateTooltip(mesh.currentGroup, params.tooltip);
        MeshBuilder.addCage(mesh, BoxState.xform, BoxCage(), params.edge_radius, 2, 8);
    }
}
