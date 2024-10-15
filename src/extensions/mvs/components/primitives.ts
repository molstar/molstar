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
import { round } from '../../../mol-util';
import { Color } from '../../../mol-util/color';
import { ParamDefinition as PD } from '../../../mol-util/param-definition';
import { capitalize } from '../../../mol-util/string';
import { rowsToExpression, rowToExpression } from '../helpers/selections';
import { collectMVSReferences, decodeColor } from '../helpers/utils';
import { ValueFor } from '../tree/generic/params-schema';
import { MVSPrimitive, MVSPrimitiveOptions, MVSPrimitiveParams } from '../tree/mvs/mvs-tree';
import { PrimitiveComponentExpressionT } from '../tree/mvs/param-types';
import { MVSTransform } from './annotation-structure-component';

type PrimitiveComponentExpression = ValueFor<typeof PrimitiveComponentExpressionT>
type MVSPositionT = [number, number, number] | PrimitiveComponentExpression | PrimitiveComponentExpression[]

export const MVSLabelProps: PD.Values<Text.Params> = {
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

export interface MVSPrimitiveBuilderContext {
    defaultStructure?: Structure;
    structureRefs: Record<string, Structure | undefined>;
    globalOptions: MVSPrimitiveOptions;
    positionCache: Map<string, Vec3>;
}

interface MeshBuilderState {
    mesh: MeshBuilder.State;
    colors: Map<number, number>;
    tooltips: Map<number, string>;
}

interface LabelBuilderState {
    group: number;
    labels: TextBuilder;
    colors: Map<number, number>;
    sizes: Map<number, number>;
}

export function getPrimitiveStructureRefs(primitives: MVSPrimitive[]) {
    const refs = new Set<string>();
    for (const p of primitives) {
        const b = Builders[p.kind];
        if (b) b[2](p, refs);
    }
    return refs;
}

function addRef(position: MVSPositionT, refs: Set<string>) {
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

function resolvePosition(context: MVSPrimitiveBuilderContext, position: MVSPositionT, target: Vec3) {
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

function buildPrimitiveMesh(context: MVSPrimitiveBuilderContext, primives: MVSPrimitive[], prev?: Mesh): Shape<Mesh> {
    const meshBuilder = MeshBuilder.createState(1024, 1024, prev);
    const state: MeshBuilderState = { mesh: meshBuilder, colors: new Map(), tooltips: new Map() };

    meshBuilder.currentGroup = -1;

    for (const p of primives) {
        const b = Builders[p.kind];
        if (!b) {
            console.warn(`Primitive ${p.kind} not supported`);
            continue;
        }
        b[0](context, state, p);
    }

    const { colors, tooltips } = state;
    const tooltip = context.globalOptions?.default_tooltip ?? '';
    const color = decodeColor(context.globalOptions?.default_color as string | undefined) ?? 0x0;

    return Shape.create(
        'Mesh',
        primives,
        MeshBuilder.getMesh(meshBuilder),
        (g) => colors.get(g) as Color ?? color as Color,
        (g) => 1,
        (g) => tooltips.get(g) ?? tooltip,
    );
}

function buildPrimitiveLabels(context: MVSPrimitiveBuilderContext, primives: MVSPrimitive[], prev?: Text): Shape<Text> {
    const labelsBuilder = TextBuilder.create(MVSLabelProps, 1024, 1024, prev);
    const state: LabelBuilderState = { group: -1, labels: labelsBuilder, colors: new Map(), sizes: new Map() };

    for (const p of primives) {
        const b = Builders[p.kind];
        if (!b) {
            console.warn(`Primitive ${p.kind} not supported`);
            continue;
        }
        b[1](context, state, p);
    }

    const { colors, sizes } = state;
    return Shape.create(
        'Labels',
        primives,
        labelsBuilder.getText(),
        (g) => colors.get(g) as Color ?? 0x0 as Color,
        (g) => sizes.get(g) ?? 1,
        (g) => '',
    );
}

const Builders: Record<MVSPrimitive['kind'], [
    mesh: (context: MVSPrimitiveBuilderContext, state: MeshBuilderState, params: any) => void,
    label: (context: MVSPrimitiveBuilderContext, state: LabelBuilderState, params: any) => void,
    resolveRefs: (params: any, refs: Set<string>) => void
]> = {
    mesh: [addMesh, noOp, noOp],
    line: [addLineMesh, noOp, resolveLineRefs],
    distance_measurement: [addDistanceMesh, addDistanceLabel, resolveLineRefs],
};

function noOp() { }

function addMesh(context: MVSPrimitiveBuilderContext, { mesh, colors, tooltips }: MeshBuilderState, params: MVSPrimitiveParams<'mesh'>) {
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

function resolveLineRefs(params: MVSPrimitiveParams<'line' | 'distance_measurement'>, refs: Set<string>) {
    addRef(params.start, refs);
    addRef(params.end, refs);
}

const lStart = Vec3.zero();
const lEnd = Vec3.zero();

function addLineMesh(context: MVSPrimitiveBuilderContext, { mesh, colors, tooltips }: MeshBuilderState, params: MVSPrimitiveParams<'line'>, options?: { skipResolvePosition?: boolean }) {
    const group = ++mesh.currentGroup;
    if (!options?.skipResolvePosition) {
        resolvePosition(context, params.start, lStart);
        resolvePosition(context, params.end, lEnd);
    }
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

function getDistanceLabel(context: MVSPrimitiveBuilderContext, params: MVSPrimitiveParams<'distance_measurement'>) {
    resolvePosition(context, params.start, lStart);
    resolvePosition(context, params.end, lEnd);

    const dist = Vec3.distance(lStart, lEnd);
    const distance = `${round(dist, 2)} Å`;
    const label = typeof params.label_template === 'string' ? params.label_template.replace('{{distance}}', distance) : distance;

    return label;
}

function addDistanceMesh(context: MVSPrimitiveBuilderContext, state: MeshBuilderState, params: MVSPrimitiveParams<'distance_measurement'>) {
    const tooltip = getDistanceLabel(context, params);
    addLineMesh(context, state, { ...params, tooltip } as any, { skipResolvePosition: true });
}

const labelPos = Vec3.zero();

function addDistanceLabel(context: MVSPrimitiveBuilderContext, state: LabelBuilderState, params: MVSPrimitiveParams<'distance_measurement'>) {
    const { labels, colors, sizes } = state;
    const group = ++state.group;
    resolvePosition(context, params.start, lStart);
    resolvePosition(context, params.end, lEnd);

    const dist = Vec3.distance(lStart, lEnd);
    const distance = `${round(dist, 2)} Å`;
    let size: number | undefined;

    if (params.label_size === 'auto') {
        size = Math.max(dist * (params.label_auto_size_scale ?? 0.2), params.label_auto_size_min ?? 0.01);
    } else if (typeof params.label_size === 'number') {
        size = params.label_size;
    }

    if (typeof size === 'number') sizes.set(group, size);
    const color = decodeColor(params.label_color);
    if (typeof color === 'number') colors.set(group, color);

    const label = typeof params.label_template === 'string' ? params.label_template.replace('{{distance}}', distance) : distance;

    Vec3.add(labelPos, lStart, lEnd);
    Vec3.scale(labelPos, labelPos, 0.5);

    labels.add(label, labelPos[0], labelPos[1], labelPos[2], 0.5, 1, group);
}



/** ========== Plugin transforms ============== */

export class MVSPrimitivesData extends SO.Create<{ primitives: MVSPrimitive[], options: MVSPrimitiveOptions, context: MVSPrimitiveBuilderContext }>({ name: 'Primitive Data', typeClass: 'Object' }) { }
export class MVSPrimitiveShapes extends SO.Create<{ mesh?: Shape<Mesh>, labels?: Shape<Text> }>({ name: 'Primitive Shapes', typeClass: 'Object' }) { }

export type MVSInlinePrimitiveData = typeof MVSInlinePrimitiveData
export const MVSInlinePrimitiveData = MVSTransform({
    name: 'mvs-inline-primitive-data',
    display: { name: 'MVS Primitives' },
    from: [SO.Root, SO.Molecule.Structure],
    to: MVSPrimitivesData,
    params: {
        primitives: PD.Value<MVSPrimitive[]>([], { isHidden: true }),
        options: PD.Value<MVSPrimitiveOptions>({} as any, { isHidden: true }),
    },
})({
    apply({ a, params }) {
        return new MVSPrimitivesData({
            primitives: params.primitives,
            options: params.options,
            context: {
                defaultStructure: SO.Molecule.Structure.is(a) ? a.data : undefined,
                structureRefs: { },
                globalOptions: params.options,
                positionCache: new Map(),
            }
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
        kind: PD.Text<'mesh' | 'labels'>('mesh')
    }
})({
    apply({ a, params, dependencies }) {
        const structureRefs = dependencies ? collectMVSReferences([SO.Molecule.Structure], dependencies) : {};
        const context: MVSPrimitiveBuilderContext = { ...a.data.context, structureRefs };

        const label = capitalize(params.kind);
        if (params.kind === 'mesh') {
            return new SO.Shape.Provider({
                label,
                data: { primitives: a.data.primitives, context },
                params: Mesh.Params,
                getShape: (_, data, __, prev: any) => buildPrimitiveMesh(data.context, data.primitives, prev),
                geometryUtils: Mesh.Utils,
            }, { label });
        } else if (params.kind === 'labels') {
            return new SO.Shape.Provider({
                label,
                data: { primitives: a.data.primitives, context },
                params: Text.Params,
                getShape: (_, data, __, prev: any) => buildPrimitiveLabels(data.context, data.primitives, prev),
                geometryUtils: Text.Utils,
            }, { label });
        }

        return StateObject.Null;
    }
});