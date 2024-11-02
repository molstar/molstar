/**
 * Copyright (c) 2024 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Lines } from '../../../mol-geo/geometry/lines/lines';
import { LinesBuilder } from '../../../mol-geo/geometry/lines/lines-builder';
import { addFixedCountDashedCylinder, addSimpleCylinder, BasicCylinderProps } from '../../../mol-geo/geometry/mesh/builder/cylinder';
import { Mesh } from '../../../mol-geo/geometry/mesh/mesh';
import { MeshBuilder } from '../../../mol-geo/geometry/mesh/mesh-builder';
import { Text } from '../../../mol-geo/geometry/text/text';
import { TextBuilder } from '../../../mol-geo/geometry/text/text-builder';
import { Box } from '../../../mol-geo/primitive/box';
import { PrimitiveBuilder } from '../../../mol-geo/primitive/primitive';
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
import { Asset } from '../../../mol-util/assets';
import { Color } from '../../../mol-util/color';
import { ParamDefinition as PD } from '../../../mol-util/param-definition';
import { capitalize } from '../../../mol-util/string';
import { rowsToExpression, rowToExpression } from '../helpers/selections';
import { collectMVSReferences, decodeColor } from '../helpers/utils';
import { MolstarNode, MolstarSubtree } from '../tree/molstar/molstar-tree';
import { MVSPrimitive, MVSPrimitiveOptions, MVSPrimitiveParams } from '../tree/mvs/mvs-primitives';
import { MVSNode } from '../tree/mvs/mvs-tree';
import { isComponentExpression, isPrimitiveComponentExpressions, isVector3, PrimitivePositionT } from '../tree/mvs/param-types';
import { MVSTransform } from './annotation-structure-component';

export function getPrimitiveStructureRefs(primitives: MolstarSubtree<'primitives'>) {
    const refs = new Set<string>();
    for (const c of primitives.children ?? []) {
        if (c.kind !== 'primitive') continue;
        const p = c.params as unknown as MVSPrimitive;
        Builders[p.kind]?.[3].refs?.(p, refs);
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
        // TODO: some list to not type everything
        kind: PD.Text<'mesh' | 'labels' | 'lines' | 'box' | 'cylinder'>('mesh')
    }
})({
    apply({ a, params, dependencies }) {
        const structureRefs = dependencies ? collectMVSReferences([SO.Molecule.Structure], dependencies) : {};
        const context: PrimitiveBuilderContext = { ...a.data, structureRefs };

        const label = capitalize(params.kind);
        // TODO: case statement

        // 
        if (params.kind === 'mesh') {
            if (!hasPrimitiveKind(a.data, 'mesh')) return StateObject.Null;

            return new SO.Shape.Provider({
                label,
                data: context,
                params: PD.withDefaults(Mesh.Params, { alpha: a.data.options?.transparency ?? 1 }),
                getShape: (_, data, __, prev: any) => buildPrimitiveMesh(data, prev?.geometry),
                geometryUtils: Mesh.Utils,
            }, { label });
        } else if (params.kind === 'labels') {
            if (!hasPrimitiveKind(a.data, 'label')) return StateObject.Null;

            return new SO.Shape.Provider({
                label,
                data: context,
                params: PD.withDefaults(DefaultLabelParams, { alpha: a.data.options?.label_transparency ?? 1 }),
                getShape: (_, data, __, prev: any) => buildPrimitiveLabels(data, prev?.geometry),
                geometryUtils: Text.Utils,
            }, { label });
        } else if (params.kind === 'lines') {
            if (!hasPrimitiveKind(a.data, 'line')) return StateObject.Null;
            debugger;
            return new SO.Shape.Provider({
                label,
                data: context,
                params: PD.withDefaults(Lines.Params, { alpha: a.data.options?.transparency ?? 1 }),
                getShape: (_, data, __, prev: any) => buildPrimitiveLines(data, prev?.geometry),
                geometryUtils: Lines.Utils,
            }, { label });
        } else if (params.kind === 'box') {
            // make it such that this is true
            // this is false, why? there are only lines
            if (!hasPrimitiveKind(a.data, 'mesh')) {
                debugger;
                return StateObject.Null;
            };
            console.log('Building primitive box');

            return new SO.Shape.Provider({
                label,
                data: context,
                params: PD.withDefaults(Mesh.Params, { alpha: a.data.options?.transparency ?? 1 }),
                getShape: (_, data, __, prev: any) => buildPrimitiveBox(data, prev?.geometry),
                // getShape: (_, data, __, prev: any) => buildPrimitiveMesh(data, prev?.geometry),
                geometryUtils: Mesh.Utils,
            }, { label });
        } else {
            throw Error(`Kind ${params.kind} is not supported yet.`);
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
    options: MVSPrimitiveOptions;
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

const Builders: Record<MVSPrimitive['kind'], [
    mesh: (context: PrimitiveBuilderContext, state: MeshBuilderState, node: MVSNode<'primitive'>, params: any) => void,
    line: (context: PrimitiveBuilderContext, state: LineBuilderState, node: MVSNode<'primitive'>, params: any) => void,
    label: (context: PrimitiveBuilderContext, state: LabelBuilderState, node: MVSNode<'primitive'>, params: any) => void,
    // box: (context: PrimitiveBuilderContext, state: BoxBuilderState, node: MVSNode<'primitive'>, params: any) => void,
    // cylinder: (context: PrimitiveBuilderContext, state: BoxBuilderState, node: MVSNode<'primitive'>, params: any) => void,
    features: {
        mesh?: boolean | ((primitive: any, context: PrimitiveBuilderContext) => boolean),
        line?: boolean | ((primitive: any, context: PrimitiveBuilderContext) => boolean),
        label?: boolean | ((primitive: any, context: PrimitiveBuilderContext) => boolean),
        // box?: boolean | ((primitive: any, context: PrimitiveBuilderContext) => boolean),
        // cylinder?: boolean | ((primitive: any, context: PrimitiveBuilderContext) => boolean),
        refs?: (params: any, refs: Set<string>) => void
    },
]> = {
    // TODO: improve impl so that the we do not have to add many noOps for the number of primitives
    mesh: [addMesh, addMeshWireframe, noOp, {
        mesh: (m: MVSPrimitiveParams<'mesh'>) => m.show_triangles ?? true,
        line: (m: MVSPrimitiveParams<'mesh'>) => m.show_wireframe ?? false,
    }],
    lines: [addMesh, addLines, noOp, { line: true }],
    line: [addLineMesh, noOp, noOp, { mesh: true, refs: resolveLineRefs }],
    label: [noOp, noOp, addPrimitiveLabel, { label: true, refs: resolveLabelRefs }],
    distance_measurement: [addDistanceMesh, noOp, addDistanceLabel, { mesh: true, label: true, refs: resolveLineRefs }],
    // TODO: probably this - try changing first noOp to addBox mesh of something
    // box: [noOp, noOp, noOp, addBox, noOp, { box: true, refs: resolveBoxRefs }],
    // cylinder: [noOp, noOp, noOp, addBox, noOp, { box: true, refs: resolveBoxRefs }],
    box: [addBoxMesh, noOp, noOp, { mesh: true, refs: resolveBoxRefs }],
    // TODO: implement
    cylinder: [noOp, noOp, noOp, { mesh: true, refs: resolveBoxRefs }]
};


function getPrimitives(primitives: MolstarSubtree<'primitives'>) {
    return (primitives.children ?? []).filter(c => c.kind === 'primitive') as unknown as MolstarNode<'primitive'>[];
}

function addRef(position: PrimitivePositionT, refs: Set<string>) {
    if (isPrimitiveComponentExpressions(position) && position.structure_ref) {
        refs.add(position.structure_ref);
    }
}

// TODO: type for this
function hasPrimitiveKind(context: PrimitiveBuilderContext, kind: 'mesh' | 'line' | 'label') {
    debugger;
    // loops over primtives in context
    // so there are no boxes in context, just lines
    for (const c of context.primitives) {
        const p = c.params as unknown as MVSPrimitive;
        const test = Builders[p.kind]?.[3]?.[kind];
        if (typeof test === 'boolean') {
            if (test) return true;
        } else if (test?.(p, context)) {
            return true;
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

function getInstances(options: MVSPrimitiveOptions | undefined): Mat4[] | undefined {
    if (!options?.instances?.length) return undefined;
    return options.instances.map(i => Mat4.fromArray(Mat4(), i, 0));
}

function buildPrimitiveMesh(context: PrimitiveBuilderContext, prev?: Mesh): Shape<Mesh> {
    const meshBuilder = MeshBuilder.createState(1024, 1024, prev);
    const state: MeshBuilderState = { groups: new GroupManager(), mesh: meshBuilder };

    meshBuilder.currentGroup = -1;

    for (const c of context.primitives) {
        const p = c.params as unknown as MVSPrimitive;
        const b = Builders[p.kind];
        if (!b) {
            console.warn(`Primitive ${p.kind} not supported`);
            continue;
        }
        b[0](context, state, c, p);
    }

    const { colors, tooltips } = state.groups;
    const tooltip = context.options?.tooltip ?? '';
    const color = decodeColor(context.options?.color) ?? 0x0;


    return Shape.create(
        'Mesh',
        {
            kind: 'mvs-primitives',
            node: context.node,
            groupToNode: state.groups.groupToNodeMap,
        },
        MeshBuilder.getMesh(meshBuilder),
        (g) => colors.get(g) as Color ?? color as Color,
        (g) => 1,
        (g) => tooltips.get(g) ?? tooltip,
        context.instances,
    );
}

function buildPrimitiveLines(context: PrimitiveBuilderContext, prev?: Lines): Shape<Lines> {
    debugger;
    const linesBuilder = LinesBuilder.create(1024, 1024, prev);
    const state: LineBuilderState = { groups: new GroupManager(), lines: linesBuilder };

    for (const c of context.primitives) {
        const p = c.params as unknown as MVSPrimitive;
        const b = Builders[p.kind];
        if (!b) {
            console.warn(`Primitive ${p.kind} not supported`);
            continue;
        }
        b[1](context, state, c, p);
    }

    const color = decodeColor(context.options?.color) ?? 0x0;
    const { colors, sizes, tooltips } = state.groups;

    return Shape.create(
        'Lines',
        {
            kind: 'mvs-primitives',
            node: context.node,
            groupToNode: state.groups.groupToNodeMap,
        },
        linesBuilder.getLines(),
        (g) => colors.get(g) as Color ?? color as Color,
        (g) => sizes.get(g) ?? 1,
        (g) => tooltips.get(g) ?? '',
        context.instances,
    );
}


// TODO: maybe use buildPrimitiveMesh directly?
function buildPrimitiveBox(context: PrimitiveBuilderContext, prev?: Mesh): Shape<Mesh> {
    debugger;
    const meshBuilder = MeshBuilder.createState(1024, 1024, prev);
    const state: MeshBuilderState = { groups: new GroupManager(), mesh: meshBuilder };

    meshBuilder.currentGroup = -1;

    for (const c of context.primitives) {
        const p = c.params as unknown as MVSPrimitive;
        const b = Builders[p.kind];
        if (!b) {
            console.warn(`Primitive ${p.kind} not supported`);
            continue;
        }
        b[0](context, state, c, p);
    }

    const { colors, tooltips } = state.groups;
    const tooltip = context.options?.tooltip ?? '';
    const color = decodeColor(context.options?.color) ?? 0x0;
    debugger;
    console.log(context, prev);
    return Shape.create(
        'Box',
        {
            kind: 'mvs-primitives',
            node: context.node,
            groupToNode: state.groups.groupToNodeMap,
        },
        MeshBuilder.getMesh(meshBuilder),
        (g) => colors.get(g) as Color ?? color as Color,
        (g) => 1,
        (g) => tooltips.get(g) ?? tooltip,
        context.instances,
    );
}


function buildPrimitiveLabels(context: PrimitiveBuilderContext, prev?: Text): Shape<Text> {
    const labelsBuilder = TextBuilder.create(BaseLabelProps, 1024, 1024, prev);
    const state: LabelBuilderState = { groups: new GroupManager(), labels: labelsBuilder };

    for (const c of context.primitives) {
        const p = c.params as unknown as MVSPrimitive;
        const b = Builders[p.kind];
        if (!b) {
            console.warn(`Primitive ${p.kind} not supported`);
            continue;
        }
        b[2](context, state, c, p);
    }

    const color = decodeColor(context.options?.label_color) ?? 0x0;
    const { colors, sizes, tooltips } = state.groups;

    return Shape.create(
        'Labels',
        {
            kind: 'mvs-primitives',
            node: context.node,
            groupToNode: state.groups.groupToNodeMap,
        },
        labelsBuilder.getText(),
        (g) => colors.get(g) as Color ?? color as Color,
        (g) => sizes.get(g) ?? 1,
        (g) => tooltips.get(g) ?? '',
        context.instances,
    );
}

function noOp() { }

function addMesh(context: PrimitiveBuilderContext, { groups, mesh }: MeshBuilderState, node: MVSNode<'primitive'>, params: MVSPrimitiveParams<'mesh'>) {
    if (!params.show_triangles) return;

    const a = Vec3.zero();
    const b = Vec3.zero();
    const c = Vec3.zero();

    const { indices, vertices, triangle_colors, triangle_groups, group_colors, group_tooltips } = params;

    const groupSet: Map<number, number> | undefined = triangle_groups?.length ? groups.allocateMany(node, triangle_groups) : undefined;

    for (let i = 0, _i = indices.length / 3; i < _i; i++) {
        if (groupSet) {
            const grp = triangle_groups![i];
            mesh.currentGroup = groupSet.get(grp)!;
            groups.updateColor(mesh.currentGroup, group_colors?.[grp] ?? params.color);
            groups.updateTooltip(mesh.currentGroup, group_tooltips?.[grp]);
        } else {
            mesh.currentGroup = groups.allocateSingle(node);
            groups.updateColor(mesh.currentGroup, triangle_colors?.[i] ?? params.color);
            groups.updateTooltip(mesh.currentGroup, params.tooltip);
        }

        Vec3.fromArray(a, vertices, 3 * indices[3 * i]);
        Vec3.fromArray(b, vertices, 3 * indices[3 * i + 1]);
        Vec3.fromArray(c, vertices, 3 * indices[3 * i + 2]);

        MeshBuilder.addTriangle(mesh, a, b, c);
    }
}

function addMeshWireframe(context: PrimitiveBuilderContext, { groups, lines }: LineBuilderState, node: MVSNode<'primitive'>, params: MVSPrimitiveParams<'mesh'>) {
    if (!params.show_wireframe) return;

    const a = Vec3.zero();
    const b = Vec3.zero();
    const c = Vec3.zero();

    const { indices, vertices, triangle_colors, triangle_groups, group_colors, group_tooltips } = params;

    const groupSet: Map<number, number> | undefined = triangle_groups?.length ? groups.allocateMany(node, triangle_groups) : undefined;
    const radius = params.wireframe_radius ?? 1;

    for (let i = 0, _i = indices.length / 3; i < _i; i++) {
        let group: number;
        if (groupSet) {
            const grp = triangle_groups![i];
            group = groupSet.get(grp)!;
            groups.updateColor(group, params.wireframe_color ?? group_colors?.[grp]);
            groups.updateTooltip(group, group_tooltips?.[grp]);
        } else {
            group = groups.allocateSingle(node);
            groups.updateColor(group, params.wireframe_color ?? triangle_colors?.[i]);
            groups.updateTooltip(group, params.tooltip);
        }

        groups.updateSize(group, radius);

        Vec3.fromArray(a, vertices, 3 * indices[3 * i]);
        Vec3.fromArray(b, vertices, 3 * indices[3 * i + 1]);
        Vec3.fromArray(c, vertices, 3 * indices[3 * i + 2]);

        lines.add(a[0], a[1], a[2], b[0], b[1], b[2], group);
        lines.add(b[0], b[1], b[2], c[0], c[1], c[2], group);
        lines.add(c[0], c[1], c[2], a[0], a[1], a[2], group);
    }
}

function addLines(context: PrimitiveBuilderContext, { groups, lines }: LineBuilderState, node: MVSNode<'primitive'>, params: MVSPrimitiveParams<'lines'>) {
    debugger;
    const a = Vec3.zero();
    const b = Vec3.zero();

    const { indices, vertices, line_colors, line_groups, group_colors, group_tooltips, group_radius } = params;

    const groupSet: Map<number, number> | undefined = line_groups?.length ? groups.allocateMany(node, line_groups) : undefined;
    const radius = params.line_radius ?? 1;

    for (let i = 0, _i = indices.length / 2; i < _i; i++) {
        let group: number;
        if (groupSet) {
            const grp = line_groups![i];
            group = groupSet.get(grp)!;
            groups.updateColor(group, group_colors?.[grp] ?? params.color);
            groups.updateTooltip(group, group_tooltips?.[grp]);
            groups.updateSize(group, group_radius?.[grp] ?? radius);
        } else {
            group = groups.allocateSingle(node);
            groups.updateColor(group, line_colors?.[i] ?? params.color);
            groups.updateSize(group, radius);
            groups.updateTooltip(group, params.tooltip);
        }

        Vec3.fromArray(a, vertices, 3 * indices[2 * i]);
        Vec3.fromArray(b, vertices, 3 * indices[2 * i + 1]);
        lines.add(a[0], a[1], a[2], b[0], b[1], b[2], group);
    }
}

function resolveLineRefs(params: MVSPrimitiveParams<'line' | 'distance_measurement'>, refs: Set<string>) {
    addRef(params.start, refs);
    addRef(params.end, refs);
}

const lStart = Vec3.zero();
const lEnd = Vec3.zero();

function addLineMesh(context: PrimitiveBuilderContext, { groups, mesh }: MeshBuilderState, node: MVSNode<'primitive'>, params: MVSPrimitiveParams<'line'>, options?: { skipResolvePosition?: boolean }) {
    if (!options?.skipResolvePosition) {
        resolveBasePosition(context, params.start, lStart);
        resolveBasePosition(context, params.end, lEnd);
    }
    const radius = params.thickness ?? 0.05;

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
        // NOTE: this: use similar
        addSimpleCylinder(mesh, lStart, lEnd, cylinderProps);
    }
}

function addBoxMesh(context: PrimitiveBuilderContext, { groups, mesh }: MeshBuilderState, node: MVSNode<'primitive'>, params: MVSPrimitiveParams<'box'>, options?: { skipResolvePosition?: boolean }) {
    // TODO: need this?
    // if (!options?.skipResolvePosition) {
    //     // TODO: need this?
    //     resolveBasePosition(context, params.start, lStart);
    //     resolveBasePosition(context, params.end, lEnd);
    // }

    const a = Vec3(), b = Vec3(), c = Vec3(), d = Vec3();
    // TODO: as_edges, rotation, edge_radius => cage, rotation matrices
    const { center, extent, box_groups, scaling, as_edges, edge_radius, rotation, translation } = params;
    const mat4 = Mat4.identity();
    const triangleCount = 12;
    // const vertexCount = perforated ? 12 * 3 : 6 * 4;
    const vertexCount = 6 * 4;
    const builder = PrimitiveBuilder(triangleCount, vertexCount);


    mesh.currentGroup = groups.allocateSingle(node);
    groups.updateColor(mesh.currentGroup, params.color);
    // groups.updateTooltip(mesh.currentGroup, params.tooltip);

    // if (params.dash_length) {
    //     const dist = Vec3.distance(lStart, lEnd);
    //     const count = Math.ceil(dist / (2 * params.dash_length));
    //     addFixedCountDashedCylinder(mesh, lStart, lEnd, 1.0, count, true, cylinderProps);
    // } else {
    //     // NOTE: this: use similar
    //     addSimpleCylinder(mesh, lStart, lEnd, cylinderProps);
    // }
    // addSimpleCylinder(mesh, lStart, lEnd, cylinderProps);

    // S1
    Vec3.set(a, 0, 0, 0);
    Vec3.set(b, 0, 0, 1);
    Vec3.set(c, 0, 1, 1);
    Vec3.set(d, 0, 1, 0);
    builder.addQuad(a, b, c, d);

    // S2
    Vec3.set(a, 0, 0, 1);
    Vec3.set(b, 0, 1, 1);
    Vec3.set(c, 1, 1, 1);
    Vec3.set(d, 1, 0, 1);
    builder.addQuad(a, b, c, d);

    // S3
    Vec3.set(a, 1, 1, 1);
    Vec3.set(b, 1, 0, 1);
    Vec3.set(c, 1, 0, 0);
    Vec3.set(d, 1, 1, 0);
    builder.addQuad(a, b, c, d);

    // S4
    Vec3.set(a, 1, 0, 0);
    Vec3.set(b, 1, 1, 0);
    Vec3.set(c, 0, 1, 0);
    Vec3.set(d, 0, 0, 0);
    builder.addQuad(a, b, c, d);


    // TOP
    Vec3.set(a, 0, 0, 0);
    Vec3.set(b, 0, 0, 1);
    Vec3.set(c, 1, 0, 1);
    Vec3.set(d, 1, 0, 0);
    builder.addQuad(a, b, c, d);

    // BOTTOM
    Vec3.set(a, 0, 1, 0);
    Vec3.set(b, 0, 1, 1);
    Vec3.set(c, 1, 1, 1);
    Vec3.set(d, 1, 1, 0);
    builder.addQuad(a, b, c, d);
    // debugger;
    // return builder.getPrimitive();
    
    const translation1 = Vec3.create(0.5, 0.5, 0.5);
    const scaling1 = Vec3.create(1, 1, 1);
    // const mat4 = Mat4.identity();
    Mat4.scale(mat4, mat4, scaling1);
    Mat4.translate(mat4, mat4, translation1);

    MeshBuilder.addPrimitive(mesh, mat4, Box());

    // const prim = builder.getPrimitive();
    // MeshBuilder.addPrimitive(prim);


    // const box = Box();
    // MeshBuilder.addPrimitive(mesh, mat4, prim);

}

function getDistanceLabel(context: PrimitiveBuilderContext, params: MVSPrimitiveParams<'distance_measurement'>) {
    resolveBasePosition(context, params.start, lStart);
    resolveBasePosition(context, params.end, lEnd);

    const dist = Vec3.distance(lStart, lEnd);
    const distance = `${round(dist, 2)} Å`;
    const label = typeof params.label_template === 'string' ? params.label_template.replace('{{distance}}', distance) : distance;

    return label;
}

function addDistanceMesh(context: PrimitiveBuilderContext, state: MeshBuilderState, node: MVSNode<'primitive'>, params: MVSPrimitiveParams<'distance_measurement'>) {
    const tooltip = getDistanceLabel(context, params);
    addLineMesh(context, state, node, { ...params, tooltip } as any, { skipResolvePosition: true });
}

const labelPos = Vec3.zero();

function addDistanceLabel(context: PrimitiveBuilderContext, state: LabelBuilderState, node: MVSNode<'primitive'>, params: MVSPrimitiveParams<'distance_measurement'>) {
    const { labels, groups } = state;
    resolveBasePosition(context, params.start, lStart);
    resolveBasePosition(context, params.end, lEnd);

    const dist = Vec3.distance(lStart, lEnd);
    const distance = `${round(dist, 2)} Å`;
    const label = typeof params.label_template === 'string' ? params.label_template.replace('{{distance}}', distance) : distance;

    let size: number | undefined;
    if (params.label_size === 'auto') {
        size = Math.max(dist * (params.label_auto_size_scale ?? 0.2), params.label_auto_size_min ?? 0.01);
    } else if (typeof params.label_size === 'number') {
        size = params.label_size;
    }

    Vec3.add(labelPos, lStart, lEnd);
    Vec3.scale(labelPos, labelPos, 0.5);

    const group = groups.allocateSingle(node);
    groups.updateColor(group, params.label_color);
    groups.updateSize(group, size);

    labels.add(label, labelPos[0], labelPos[1], labelPos[2], 1.05 * (params.thickness ?? 0.05), 1, group);
}

function resolveLabelRefs(params: MVSPrimitiveParams<'label'>, refs: Set<string>) {
    addRef(params.position, refs);
}

// TODO: type for Vec3 float list
function resolveBoxRefs(params: MVSPrimitiveParams<'box'>, refs: Set<string>) {
    addRef(params.center as [number, number, number], refs);
}

function addPrimitiveLabel(context: PrimitiveBuilderContext, state: LabelBuilderState, node: MVSNode<'primitive'>, params: MVSPrimitiveParams<'label'>) {
    const { labels, groups } = state;
    resolveBasePosition(context, params.position, labelPos);

    const group = groups.allocateSingle(node);
    groups.updateColor(group, params.label_color);
    groups.updateSize(group, params.label_size);

    labels.add(params.text, labelPos[0], labelPos[1], labelPos[2], params.label_offset ?? 0, 1, group);
}
