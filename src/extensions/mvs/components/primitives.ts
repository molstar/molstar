/**
 * Copyright (c) 2024 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { addFixedCountDashedCylinder, addSimpleCylinder, BasicCylinderProps } from '../../../mol-geo/geometry/mesh/builder/cylinder';
import { Mesh } from '../../../mol-geo/geometry/mesh/mesh';
import { MeshBuilder } from '../../../mol-geo/geometry/mesh/mesh-builder';
import { Text } from '../../../mol-geo/geometry/text/text';
import { TextBuilder } from '../../../mol-geo/geometry/text/text-builder';
import { Vec3 } from '../../../mol-math/linear-algebra';
import { Loci } from '../../../mol-model/loci';
import { Shape } from '../../../mol-model/shape';
import { Structure, StructureSelection } from '../../../mol-model/structure';
import { StructureQueryHelper } from '../../../mol-plugin-state/helpers/structure-query';
import { PluginStateObject as SO } from '../../../mol-plugin-state/objects';
import { Expression } from '../../../mol-script/language/expression';
import { StateObject } from '../../../mol-state';
import { Color } from '../../../mol-util/color';
import { ParamDefinition as PD } from '../../../mol-util/param-definition';
import { capitalize } from '../../../mol-util/string';
import { rowsToExpression, rowToExpression } from '../helpers/selections';
import { collectMVSReferences, decodeColor } from '../helpers/utils';
import { ValueFor } from '../tree/generic/params-schema';
import { MVSNode } from '../tree/mvs/mvs-tree';
import { ComponentExpressionT } from '../tree/mvs/param-types';
import { MVSTransform } from './annotation-structure-component';

export interface PrimitiveComponentExpression extends ValueFor<typeof ComponentExpressionT> {
    structure_ref?: string;
}

export type PositionT = [number, number, number] | PrimitiveComponentExpression | PrimitiveComponentExpression[]

export interface PrimitiveOptions {
    color?: string;
    transparency?: number;
    tooltip?: string;
}

export type Primitive =
    | { kind: 'primitive_mesh', params: MVSNode<'primitive_mesh'>['params'] }
    | { kind: 'primitive_line', params: MVSNode<'primitive_line'>['params'] }
    | { kind: 'primitive_distance_measurement', params: MVSNode<'primitive_distance_measurement'>['params'] };

export interface PrimitiveBuilderContext {
    mesh?: Mesh;
    labels?: Text;

    defaultStructure?: Structure;
    structureRefs: Record<string, Structure | undefined>;
    globalOptions: PrimitiveOptions;
    positionCache: Map<string, Vec3>;
}

interface BuilderState {
    mesh: MeshBuilder.State;
    labels: TextBuilder;
    colors: Map<number, number>;
    tooltips: Map<number, string>;
}

export function getPrimitiveStructureRefs(primitives: Primitive[]) {
    const refs = new Set<string>();
    for (const p of primitives) {
        const b = Builders[p.kind];
        if (b) b[1](p, refs);
    }
    return refs;
}

function addRef(position: PositionT, refs: Set<string>) {
    if (Array.isArray(position)) {
        if (typeof position[0] === 'number') {
            return;
        }
        for (const p of position as PrimitiveComponentExpression[]) {
            if (p.structure_ref) refs.add(p.structure_ref);
        }
    }
    if ((position as PrimitiveComponentExpression).structure_ref) {
        refs.add((position as PrimitiveComponentExpression).structure_ref!);
    }
}

function resolvePosition(context: PrimitiveBuilderContext, position: PositionT, target: Vec3) {
    let expr: Expression;
    let pivotRef: string | undefined;
    if (Array.isArray(position)) {
        if (typeof position[0] === 'number') {
            return Vec3.copy(target, position as Vec3);
        }
        if (position.length === 0) return Vec3.set(target, 0, 0, 0);

        const exprs = position as PrimitiveComponentExpression[];
        pivotRef = exprs[0].structure_ref;
        for (const e of exprs) {
            if (pivotRef !== e.structure_ref) throw new Error('All position expressions must point to the same structure');
        }

        expr = rowsToExpression(position as any);
    } else {
        expr = rowToExpression(position as any);
    }

    const pivot = !pivotRef ? context.defaultStructure : context.structureRefs[pivotRef];
    if (!pivot) {
        throw new Error(`Structure with ref '${pivotRef ?? '<default>'}' not found.`);
    }

    if (!context.defaultStructure) return Vec3.set(target, 0, 0, 0);

    const cackeKey = JSON.stringify(position);
    if (context.positionCache.has(cackeKey)) {
        return Vec3.copy(target, context.positionCache.get(cackeKey)!);
    }

    const { selection } = StructureQueryHelper.createAndRun(context.defaultStructure, expr);

    if (StructureSelection.isEmpty(selection)) {
        Vec3.set(target, 0, 0, 0);
    } else {
        const loci = StructureSelection.toLociWithSourceUnits(selection);
        Loci.getCenter(loci, target);
    }

    context.positionCache.set(cackeKey, Vec3.clone(target));
}

const LabelProps: PD.Values<Text.Params> = {
    ...PD.getDefaultValues(Text.Params),
    attachment: 'middle-center',
    fontQuality: 3,
    fontWeight: 'normal',
    borderWidth: 0.3,
    borderColor: Color(0x0),
    background: false,
    backgroundOpacity: 0.5,
    tether: false,
};

function buildPrimitiveShapes(context: PrimitiveBuilderContext, primives: Primitive[]): { mesh?: Shape<Mesh>, labels?: Shape<Text> } {
    const meshBuilder = MeshBuilder.createState(1024, 1024, context.mesh);
    const labelsBuilder = TextBuilder.create(LabelProps, 1024, 1024, context.labels);
    const state: BuilderState = { mesh: meshBuilder, labels: labelsBuilder, colors: new Map(), tooltips: new Map() };

    meshBuilder.currentGroup = -1;

    for (const p of primives) {
        const b = Builders[p.kind];
        if (!b) throw new Error(`Primitive ${p.kind} not supported`);
        b[0](context, state, p.params as any);
    }

    const { colors, tooltips } = state;
    const tooltip = context.globalOptions.tooltip ?? '';
    const color = decodeColor(context.globalOptions.color) ?? 0x0;

    const mesh = meshBuilder.currentGroup !== -1 ? Shape.create(
        'Mesh',
        primives,
        MeshBuilder.getMesh(meshBuilder),
        (g) => colors.get(g) as Color ?? color as Color,
        () => 1,
        (g) => tooltips.get(g) ?? tooltip,
    ) : undefined;

    // TODO
    const labels = !labelsBuilder.isEmpty ? Shape.create(
        'Labels',
        primives,
        labelsBuilder.getText(),
        (g) => color as Color,
        () => 1,
        (g) => '',
    ) : undefined;

    return { mesh, labels };
}

const Builders: Record<Primitive['kind'], [
    build: (context: PrimitiveBuilderContext, state: BuilderState, params: any) => void,
    resolveRefs: (params: any, refs: Set<string>) => void
]> = {
    primitive_mesh: [addMesh, () => {}],
    primitive_line: [addLine, (params: MVSNode<'primitive_line'>['params'], refs) => {
        addRef(params.start, refs);
        addRef(params.end, refs);
    }],
    primitive_distance_measurement: [addLine, (params: MVSNode<'primitive_distance_measurement'>['params'], refs) => {
        addRef(params.start, refs);
        addRef(params.end, refs);
    }],
};

function addMesh(context: PrimitiveBuilderContext, { mesh, colors, tooltips }: BuilderState, params: MVSNode<'primitive_mesh'>['params']) {
    const a = Vec3.zero();
    const b = Vec3.zero();
    const c = Vec3.zero();

    const { indices, vertices, triangle_colors, triangle_groups, group_colors, group_tooltips } = params;

    const baseGroup = mesh.currentGroup;
    const groupOffsets = new Map<number, number>();
    if (triangle_groups) {
        for (const g of triangle_groups) {
            if (groupOffsets.has(g)) continue;
            groupOffsets.set(g, groupOffsets.size + 1);
        }
        mesh.currentGroup += groupOffsets.size + 1;
    }

    for (let i = 0, _i = indices.length / 3; i < _i; i++) {
        let group: number;
        let color: number | undefined;
        let tooltip: string | undefined;

        if (triangle_groups) {
            const grp = triangle_groups[i];
            mesh.currentGroup = baseGroup + groupOffsets.get(grp)!;
            group = mesh.currentGroup;
            color = decodeColor(group_colors?.[grp]);
            tooltip = group_tooltips?.[grp];
        } else {
            group = ++mesh.currentGroup;
            color = decodeColor(triangle_colors?.[i]);
        }

        if (typeof color !== 'undefined') colors.set(group, color);
        if (tooltip) tooltips.set(group, tooltip);

        Vec3.fromArray(a, vertices, 3 * indices[3 * i]);
        Vec3.fromArray(b, vertices, 3 * indices[3 * i + 1]);
        Vec3.fromArray(c, vertices, 3 * indices[3 * i + 2]);

        MeshBuilder.addTriangle(mesh, a, b, c);
    }
}

const lStart = Vec3.zero();
const lEnd = Vec3.zero();

function addLine(context: PrimitiveBuilderContext, { mesh, colors, tooltips }: BuilderState, params: MVSNode<'primitive_line'>['params']) {
    const group = ++mesh.currentGroup;
    resolvePosition(context, params.start, lStart);
    resolvePosition(context, params.end, lEnd);
    const radius = params.thickness ?? 0.05;

    const cylinderProps: BasicCylinderProps = {
        radiusBottom: radius,
        radiusTop: radius,
        topCap: true,
        bottomCap: true,
    };

    if (params.dash_length) {
        // TODO: support other dash params
        const dist = Vec3.distance(lStart, lEnd);
        const count = Math.ceil(dist / (2 * params.dash_length));
        addFixedCountDashedCylinder(mesh, lStart, lEnd, 1.0, count, true, cylinderProps);
    } else {
        addSimpleCylinder(mesh, lStart, lEnd, cylinderProps);
    }

    const color = decodeColor(params?.color as string);
    if (typeof color !== 'undefined') colors.set(group, color);
    if (params.tooltip) tooltips.set(group, params.tooltip);
}

/** ========== Plugin transforms ============== */

export class MVSPrimitivesData extends SO.Create<{ structure?: Structure, primitives: Primitive[], options: PrimitiveOptions }>({ name: 'Primitive Data', typeClass: 'Object' }) { }
export class MVSPrimitiveShapes extends SO.Create<{ mesh?: Shape<Mesh>, labels?: Shape<Text> }>({ name: 'Primitive Shapes', typeClass: 'Object' }) { }

export type MVSInlinePrimitiveData = typeof MVSInlinePrimitiveData
export const MVSInlinePrimitiveData = MVSTransform({
    name: 'mvs-inline-primitive-data',
    display: { name: 'MVS Primitives' },
    from: [SO.Root, SO.Molecule.Structure],
    to: MVSPrimitivesData,
    params: {
        primitives: PD.Value<Primitive[]>([], { isHidden: true }),
        options: PD.Value<PrimitiveOptions>({} as any, { isHidden: true }),
    },
})({
    apply({ a, params }) {
        return new MVSPrimitivesData({
            structure: SO.Molecule.Structure.is(a) ? a.data : undefined,
            primitives: params.primitives,
            options: params.options,
        }, { label: 'Primitive Data' });
    }
})

export type MVSBuildPrimitiveShapes = typeof MVSBuildPrimitiveShapes
export const MVSBuildPrimitiveShapes = MVSTransform({
    name: 'mvs-build-primitive-shapes',
    display: { name: 'MVS Primitives' },
    from: MVSPrimitivesData,
    to: MVSPrimitiveShapes,
})({
    apply({ a, dependencies }) {
        const structureRefs = dependencies ? collectMVSReferences([SO.Molecule.Structure], dependencies) : {};
        const context: PrimitiveBuilderContext = {
            defaultStructure: a.data.structure,
            structureRefs,
            globalOptions: a.data.options,
            positionCache: new Map(),
        };
        const shapes = buildPrimitiveShapes(context, a.data.primitives);
        return new MVSPrimitiveShapes(shapes, { label: 'Primitive Shapes ' });
    }
});

export type MVSBuildPrimitiveShapeProvider = typeof MVSBuildPrimitiveShapes
export const MVSBuildPrimitiveShapeProvider = MVSTransform({
    name: 'mvs-build-primitive-shape-provider',
    display: { name: 'MVS Primitives' },
    from: MVSPrimitiveShapes,
    to: SO.Shape.Provider,
    params: {
        kind: PD.Text<'mesh' | 'labels'>('mesh')
    }
})({
    apply({ a, params }) {
        if (!a.data[params.kind]) return StateObject.Null;
        const geo = params.kind === 'mesh' ? Mesh : Text;
        return new SO.Shape.Provider({
            label: capitalize(params.kind),
            data: a.data[params.kind],
            params: geo.Params,
            getShape: (_, data: any) => data,
            geometryUtils: geo.Utils
        }, { label: capitalize(params.kind) });
    }
});