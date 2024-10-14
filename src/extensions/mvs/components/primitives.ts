/**
 * Copyright (c) 2024 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { addSimpleCylinder } from '../../../mol-geo/geometry/mesh/builder/cylinder';
import { Mesh } from '../../../mol-geo/geometry/mesh/mesh';
import { MeshBuilder } from '../../../mol-geo/geometry/mesh/mesh-builder';
import { Vec3 } from '../../../mol-math/linear-algebra';
import { Shape } from '../../../mol-model/shape';
import { Structure } from '../../../mol-model/structure';
import { PluginStateObject as SO } from '../../../mol-plugin-state/objects';
import { Color } from '../../../mol-util/color';
import { ParamDefinition as PD } from '../../../mol-util/param-definition';
import { decodeColor } from '../helpers/utils';
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
    | { kind: 'mesh', params: MVSNode<'mesh'>['params'] }
    | { kind: 'line', params: MVSNode<'line'>['params'] };

export interface PrimitiveBuilderContext {
    mesh?: Mesh;
    defaultStructure?: Structure;
    structureRefs: Record<string, Structure | undefined>;
    globalOptions: PrimitiveOptions;
}

interface BuilderState {
    builder: MeshBuilder.State;
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
    if (Array.isArray(position)) {
        if (typeof position[0] === 'number') {
            return Vec3.copy(target, position as Vec3);
        }
    }

    throw new Error('not yet implemented');
}

function buildPrimitiveMesh(context: PrimitiveBuilderContext, primives: Primitive[]): Shape<Mesh> {
    const builder = MeshBuilder.createState(1024, 1024, context.mesh);
    const state: BuilderState = { builder, colors: new Map(), tooltips: new Map() };
    builder.currentGroup = -1;

    for (const p of primives) {
        const b = Builders[p.kind];
        if (!b) throw new Error(`Primitive ${p.kind} not supported`);
        b[0](context, state, p.params as any);
    }

    const { colors, tooltips } = state;
    const tooltip = context.globalOptions.tooltip ?? '';
    const color = decodeColor(context.globalOptions.color) ?? 0x0;

    return Shape.create(
        'Mesh',
        primives,
        MeshBuilder.getMesh(builder),
        (g) => colors.get(g) as Color ?? color as Color,
        () => 1,
        (g) => tooltips.get(g) ?? tooltip,
    );
}

const Builders: Record<Primitive['kind'], [
    build: (context: PrimitiveBuilderContext, state: BuilderState, params: any) => void,
    resolveRefs: (params: any, refs: Set<string>) => void
]> = {
    mesh: [addMesh, () => {}],
    line: [addLine, (params: MVSNode<'line'>['params'], refs) => {
        addRef(params.start, refs);
        addRef(params.end, refs);
    }],
};

function addMesh(context: PrimitiveBuilderContext, { builder, colors, tooltips }: BuilderState, params: MVSNode<'mesh'>['params']) {
    const a = Vec3.zero();
    const b = Vec3.zero();
    const c = Vec3.zero();

    const { indices, vertices, triangle_colors, triangle_groups, group_colors, group_tooltips } = params;

    const baseGroup = builder.currentGroup;
    const groupOffsets = new Map<number, number>();
    if (triangle_groups) {
        for (const g of triangle_groups) {
            if (groupOffsets.has(g)) continue;
            groupOffsets.set(g, groupOffsets.size + 1);
        }
        builder.currentGroup += groupOffsets.size + 1;
    }


    for (let i = 0, _i = indices.length / 3; i < _i; i++) {
        let group: number;
        let color: number | undefined;
        let tooltip: string | undefined;

        if (triangle_groups) {
            const grp = triangle_groups[i];
            builder.currentGroup = baseGroup + groupOffsets.get(grp)!;
            group = builder.currentGroup;
            color = decodeColor(group_colors?.[grp]);
            tooltip = group_tooltips?.[grp];
        } else {
            group = ++builder.currentGroup;
            color = decodeColor(triangle_colors?.[i]);
        }

        if (typeof color !== 'undefined') colors.set(group, color);
        if (tooltip) tooltips.set(group, tooltip);

        Vec3.fromArray(a, vertices, 3 * indices[3 * i]);
        Vec3.fromArray(b, vertices, 3 * indices[3 * i + 1]);
        Vec3.fromArray(c, vertices, 3 * indices[3 * i + 2]);

        MeshBuilder.addTriangle(builder, a, b, c);
    }
}

const lStart = Vec3.zero();
const lEnd = Vec3.zero();

function addLine(context: PrimitiveBuilderContext, { builder, colors, tooltips }: BuilderState, params: MVSNode<'line'>['params']) {
    const group = ++builder.currentGroup;
    resolvePosition(context, params.start, lStart);
    resolvePosition(context, params.end, lEnd);
    const radius = params.thickness ?? 0.05;
    // TODO: support dashes etc
    addSimpleCylinder(builder, lStart, lEnd, { radiusBottom: radius, radiusTop: radius, topCap: true, bottomCap: true });
    const color = decodeColor(params?.color as string);
    if (typeof color !== 'undefined') colors.set(group, color);
    if (params.tooltip) tooltips.set(group, params.tooltip);
}

export type MVSInlinePrimitives = typeof MVSInlinePrimitives
export const MVSInlinePrimitives = MVSTransform({ // TODO: move MVSTransform to a common file
    name: 'mvs-inline-primitives',
    display: { name: 'MVS Primitives' },
    from: [SO.Root, SO.Molecule.Structure],
    to: SO.Shape.Provider,
    params: {
        primitives: PD.Value<Primitive[]>([], { isHidden: true }),
        options: PD.Value<PrimitiveOptions>({} as any, { isHidden: true }),
    },
})({
    apply({ params, dependencies }) {
        return new SO.Shape.Provider({
            label: 'Primitives',
            data: params,
            params: Mesh.Params,
            getShape: (_, data: typeof params) => {
                const mesh = buildPrimitiveMesh({
                    structureRefs: {},
                    globalOptions: data.options,
                }, data.primitives);
                return mesh;
            },
            geometryUtils: Mesh.Utils
        }, { label: 'Box' });
    }
});