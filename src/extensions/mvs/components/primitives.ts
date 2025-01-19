/**
 * Copyright (c) 2024 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Adam Midlik <midlik@gmail.com>
 */

import { Lines } from '../../../mol-geo/geometry/lines/lines';
import { LinesBuilder } from '../../../mol-geo/geometry/lines/lines-builder';
import { addFixedCountDashedCylinder, addSimpleCylinder, BasicCylinderProps } from '../../../mol-geo/geometry/mesh/builder/cylinder';
import { Mesh } from '../../../mol-geo/geometry/mesh/mesh';
import { MeshBuilder } from '../../../mol-geo/geometry/mesh/mesh-builder';
import { Text } from '../../../mol-geo/geometry/text/text';
import { TextBuilder } from '../../../mol-geo/geometry/text/text-builder';
import { Box3D, Sphere3D } from '../../../mol-math/geometry';
import { Mat4, Vec3 } from '../../../mol-math/linear-algebra';
import { Shape } from '../../../mol-model/shape';
import { Structure, StructureElement, StructureSelection } from '../../../mol-model/structure';
import { StructureQueryHelper } from '../../../mol-plugin-state/helpers/structure-query';
import { PluginStateObject as SO } from '../../../mol-plugin-state/objects';
import { PluginContext } from '../../../mol-plugin/context';
import { Expression } from '../../../mol-script/language/expression';
import { StateObject } from '../../../mol-state';
import { Task } from '../../../mol-task';
import { round } from '../../../mol-util';
import { range } from '../../../mol-util/array';
import { Asset } from '../../../mol-util/assets';
import { Color } from '../../../mol-util/color';
import { ParamDefinition as PD } from '../../../mol-util/param-definition';
import { capitalize } from '../../../mol-util/string';
import { rowsToExpression, rowToExpression } from '../helpers/selections';
import { collectMVSReferences, decodeColor } from '../helpers/utils';
import { MolstarNode, MolstarSubtree } from '../tree/molstar/molstar-tree';
import { MVSNode } from '../tree/mvs/mvs-tree';
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
            const node = JSON.parse(asset.data) as MolstarSubtree<'primitives'>;
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

export type MVSInlinePrimitiveData = typeof MVSInlinePrimitiveData
export const MVSInlinePrimitiveData = MVSTransform({
    name: 'mvs-inline-primitive-data',
    display: { name: 'MVS Primitives' },
    from: [SO.Root, SO.Molecule.Structure],
    to: MVSPrimitivesData,
    params: {
        node: PD.Value<MolstarSubtree<'primitives'>>(undefined as any, { isHidden: true }),
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

        const label = capitalize(params.kind);
        if (params.kind === 'mesh') {
            if (!hasPrimitiveKind(a.data, 'mesh')) return StateObject.Null;

            return new SO.Shape.Provider({
                label,
                data: context,
                params: PD.withDefaults(Mesh.Params, { alpha: a.data.options?.opacity ?? 1 }),
                getShape: (_, data, __, prev: any) => buildPrimitiveMesh(data, prev?.geometry),
                geometryUtils: Mesh.Utils,
            }, { label });
        } else if (params.kind === 'labels') {
            if (!hasPrimitiveKind(a.data, 'label')) return StateObject.Null;

            return new SO.Shape.Provider({
                label,
                data: context,
                params: PD.withDefaults(DefaultLabelParams, { alpha: a.data.options?.label_opacity ?? 1 }),
                getShape: (_, data, __, prev: any) => buildPrimitiveLabels(data, prev?.geometry),
                geometryUtils: Text.Utils,
            }, { label });
        } else if (params.kind === 'lines') {
            if (!hasPrimitiveKind(a.data, 'line')) return StateObject.Null;

            return new SO.Shape.Provider({
                label,
                data: context,
                params: PD.withDefaults(Lines.Params, { alpha: a.data.options?.opacity ?? 1 }),
                getShape: (_, data, __, prev: any) => buildPrimitiveLines(data, prev?.geometry),
                geometryUtils: Lines.Utils,
            }, { label });
        }

        return StateObject.Null;
    }
});

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
    positionCache: Map<string, [Sphere3D, Box3D]>;
    instances: Mat4[] | undefined;
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

function resolveBasePosition(context: PrimitiveBuilderContext, position: PrimitivePositionT, targetPosition: Vec3) {
    return resolvePosition(context, position, targetPosition, undefined, undefined);
}

const _EmptySphere = Sphere3D.zero();
const _EmptyBox = Box3D.zero();

function resolvePosition(context: PrimitiveBuilderContext, position: PrimitivePositionT, targetPosition: Vec3 | undefined, targetSphere: Sphere3D | undefined, targetBox: Box3D | undefined) {
    let expr: Expression | undefined;
    let pivotRef: string | undefined;

    if (isVector3(position)) {
        if (targetPosition) Vec3.copy(targetPosition, position as any);
        if (targetSphere) Sphere3D.set(targetSphere, position as any, 0);
        if (targetBox) Box3D.set(targetBox, position as any, position as any);
        return;
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

    const cackeKey = JSON.stringify(position);
    if (context.positionCache.has(cackeKey)) {
        const cached = context.positionCache.get(cackeKey)!;
        if (targetPosition) Vec3.copy(targetPosition, cached[0].center);
        if (targetSphere) Sphere3D.copy(targetSphere, cached[0]);
        if (targetBox) Box3D.copy(targetBox, cached[1]);
        return;
    }

    const { selection } = StructureQueryHelper.createAndRun(pivot, expr);

    let box: Box3D;
    let sphere: Sphere3D;

    if (StructureSelection.isEmpty(selection)) {
        if (targetPosition) Vec3.set(targetPosition, 0, 0, 0);
        box = _EmptyBox;
        sphere = _EmptySphere;
    } else {
        const loci = StructureSelection.toLociWithSourceUnits(selection);
        const boundary = StructureElement.Loci.getBoundary(loci);
        if (targetPosition) Vec3.copy(targetPosition, boundary.sphere.center);
        box = boundary.box;
        sphere = boundary.sphere;
    }

    if (targetSphere) Sphere3D.copy(targetSphere, sphere);
    if (targetBox) Box3D.copy(targetBox, box);

    context.positionCache.set(cackeKey, [sphere, box]);
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

function buildPrimitiveLabels(context: PrimitiveBuilderContext, prev?: Text): Shape<Text> {
    const labelsBuilder = TextBuilder.create(BaseLabelProps, 1024, 1024, prev);
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
        resolveBasePosition(context, params.start, lStart);
        resolveBasePosition(context, params.end, lEnd);
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

function getDistanceLabel(context: PrimitiveBuilderContext, params: PrimitiveParams<'distance_measurement'>) {
    resolveBasePosition(context, params.start, lStart);
    resolveBasePosition(context, params.end, lEnd);

    const dist = Vec3.distance(lStart, lEnd);
    const distance = `${round(dist, 2)} Å`;
    const label = typeof params.label_template === 'string' ? params.label_template.replace('{{distance}}', distance) : distance;

    return label;
}

function addDistanceMesh(context: PrimitiveBuilderContext, state: MeshBuilderState, node: MVSNode<'primitive'>, params: PrimitiveParams<'distance_measurement'>) {
    const tooltip = getDistanceLabel(context, params);
    addTubeMesh(context, state, node, { ...params, tooltip } as any, { skipResolvePosition: true });
}

const labelPos = Vec3.zero();

function addDistanceLabel(context: PrimitiveBuilderContext, state: LabelBuilderState, node: MVSNode<'primitive'>, params: PrimitiveParams<'distance_measurement'>) {
    const { labels, groups } = state;
    resolveBasePosition(context, params.start, lStart);
    resolveBasePosition(context, params.end, lEnd);

    const dist = Vec3.distance(lStart, lEnd);
    const distance = `${round(dist, 2)} Å`;
    const label = typeof params.label_template === 'string' ? params.label_template.replace('{{distance}}', distance) : distance;

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

    labels.add(label, labelPos[0], labelPos[1], labelPos[2], 1.05 * (params.radius), 1, group);
}

function resolveLabelRefs(params: PrimitiveParams<'label'>, refs: Set<string>) {
    addRef(params.position, refs);
}

function addPrimitiveLabel(context: PrimitiveBuilderContext, state: LabelBuilderState, node: MVSNode<'primitive'>, params: PrimitiveParams<'label'>) {
    const { labels, groups } = state;
    resolveBasePosition(context, params.position, labelPos);

    const group = groups.allocateSingle(node);
    groups.updateColor(group, params.label_color);
    groups.updateSize(group, params.label_size);

    labels.add(params.text, labelPos[0], labelPos[1], labelPos[2], params.label_offset, 1, group);
}
