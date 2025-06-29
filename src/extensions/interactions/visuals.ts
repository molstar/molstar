/**
 * Copyright (c) 2025 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { addFixedCountDashedCylinder, addSimpleCylinder, BasicCylinderProps } from '../../mol-geo/geometry/mesh/builder/cylinder';
import { Mesh } from '../../mol-geo/geometry/mesh/mesh';
import { MeshBuilder } from '../../mol-geo/geometry/mesh/mesh-builder';
import { Sphere3D } from '../../mol-math/geometry';
import { Vec3 } from '../../mol-math/linear-algebra';
import { Shape } from '../../mol-model/shape';
import { StructureElement } from '../../mol-model/structure';
import { addLinkCylinderMesh, AddLinkOptions, AddLinkParams, DefaultLinkCylinderProps, LinkStyle } from '../../mol-repr/structure/visual/util/link';
import { Color } from '../../mol-util/color';
import { ParamDefinition as PD } from '../../mol-util/param-definition';
import { stringToWords } from '../../mol-util/string';
import { InteractionKinds, StructureInteractions } from './model';

function visualParams({ color, style = 'dashed', radius = 0.04 }: { color: Color, style?: 'dashed' | 'solid', radius?: number }) {
    return PD.Group({
        color: PD.Color(color),
        style: PD.Select<'dashed' | 'solid'>(style, [['dashed', 'Dashed'], ['solid', 'Solid']]),
        radius: PD.Numeric(radius, { min: 0.01, max: 1, step: 0.01 }),
    });
}

function hydrogenVisualParams({ color, style = 'dashed', radius = 0.04, showArrow = true, arrowOffset = 0.18 }: { color: Color, style?: 'dashed' | 'solid', radius?: number, showArrow?: boolean, arrowOffset?: number }) {
    return PD.Group({
        color: PD.Color(color),
        style: PD.Select<'dashed' | 'solid'>(style, [['dashed', 'Dashed'], ['solid', 'Solid']]),
        radius: PD.Numeric(radius, { min: 0.01, max: 1, step: 0.01 }),
        showArrow: PD.Boolean(showArrow),
        arrowOffset: PD.Numeric(arrowOffset, { min: 0, max: 1, step: 0.001 }),
    });
}

export const InteractionVisualParams = {
    kinds: PD.MultiSelect(InteractionKinds, InteractionKinds.map(k => [k, stringToWords(k)])),
    styles: PD.Group({
        'unknown': visualParams({ color: Color(0x0) }),
        'ionic': visualParams({ color: Color(0xADD8E6) }),
        'pi-stacking': visualParams({ color: Color(0x1E3F66) }),
        'cation-pi': visualParams({ color: Color(0x06402B) }),
        'halogen-bond': visualParams({ color: Color(0xFFDE21) }),
        'hydrogen-bond': hydrogenVisualParams({ color: Color(0x0), style: 'solid' }),
        'weak-hydrogen-bond': hydrogenVisualParams({ color: Color(0x0) }),
        'hydrophobic': visualParams({ color: Color(0x555555) }),
        'metal-coordination': visualParams({ color: Color(0x952e8f) }),
        'covalent': PD.Group({
            color: PD.Color(Color(0x999999)),
            radius: PD.Numeric(0.1, { min: 0.01, max: 1, step: 0.01 }),
        }),
    })
};

export type InteractionVisualParams = PD.Values<typeof InteractionVisualParams>;

export function buildInteractionsShape(interactions: StructureInteractions, params: InteractionVisualParams, prev?: Mesh): Shape<Mesh> {
    const mesh = MeshBuilder.createState(interactions.elements.length * 128, 1024, prev);

    mesh.currentGroup = -1;
    const tooltips = new Map<number, string>();

    const visible = new Set(params.kinds);
    const kindsToWords = new Map(InteractionKinds.map(k => [k, stringToWords(k)]));

    const colors = new Map<number, Color>();

    const bA = { sphere: Sphere3D.zero() };
    const bB = { sphere: Sphere3D.zero() };

    const pA = Vec3();
    const pB = Vec3();

    const dir = Vec3();
    const capPos = Vec3();

    const addLinkOptions: AddLinkOptions = {
        builderState: mesh,
        props: { ...DefaultLinkCylinderProps },
    };
    const addLinkParams: AddLinkParams = {
        a: pA,
        b: pB,
        group: 0,
        linkStub: false,
        linkStyle: LinkStyle.Solid,
        linkRadius: 0,
    };

    for (const interaction of interactions.elements) {
        mesh.currentGroup++;
        if (!visible.has(interaction.info.kind)) continue;

        let tooltip: string;
        if (interaction.info.kind === 'covalent') {
            if (interaction.info.degree === 'aromatic') tooltip = 'Aromatic';
            else if (interaction.info.degree === 1) tooltip = 'Single';
            else if (interaction.info.degree === 2) tooltip = 'Double';
            else if (interaction.info.degree === 3) tooltip = 'Triple';
            else if (interaction.info.degree === 4) tooltip = 'Quadruple';
            else tooltip = 'Covalent';
        } else {
            tooltip = kindsToWords.get(interaction.info.kind) ?? interaction.info.kind;
        }
        if (interaction.sourceSchema?.description) {
            tooltip += ` (${interaction.sourceSchema.description})`;
        }
        tooltips.set(mesh.currentGroup, tooltip);

        const options = params.styles[interaction.info.kind];
        let style: 'dashed' | 'solid' = 'solid';
        if (interaction.info.kind !== 'covalent') {
            style = params.styles[interaction.info.kind].style;
        }

        colors.set(mesh.currentGroup, params.styles[interaction.info.kind].color);

        StructureElement.Loci.getBoundary(interaction.a, undefined, bA);
        StructureElement.Loci.getBoundary(interaction.b, undefined, bB);

        Vec3.sub(dir, bB.sphere.center, bA.sphere.center);
        Vec3.normalize(dir, dir);

        Vec3.copy(pA, bA.sphere.center);
        Vec3.copy(pB, bB.sphere.center);

        if (interaction.info.kind === 'hydrogen-bond' || interaction.info.kind === 'weak-hydrogen-bond') {
            const hydrogenStyle = params.styles[interaction.info.kind];
            if (hydrogenStyle.showArrow && hydrogenStyle.arrowOffset > 0) {
                Vec3.scaleAndAdd(pB, pB, dir, -hydrogenStyle.arrowOffset);
            }

            if (hydrogenStyle.showArrow) {
                const height = options.radius * 3;
                Vec3.scaleAndAdd(capPos, pB, dir, -height);
                cylinder(mesh, pA, capPos, options.radius, style);
                addSimpleCylinder(
                    mesh,
                    capPos,
                    pB,
                    { radiusTop: 0, radiusBottom: height, topCap: false, bottomCap: true }
                );
            } else {
                cylinder(mesh, pA, pB, options.radius, style);
            }
        } else {
            if (interaction.info.kind !== 'covalent') {
                cylinder(mesh, pA, pB, options.radius, style);
            } else {
                addLinkParams.group = mesh.currentGroup;
                addLinkParams.linkRadius = options.radius;
                const degree = interaction.info.degree ?? 1;
                if (degree === 'aromatic') addLinkParams.linkStyle = LinkStyle.Aromatic;
                else if (degree === 2) addLinkParams.linkStyle = LinkStyle.Double;
                else if (degree === 3) addLinkParams.linkStyle = LinkStyle.Triple;
                else addLinkParams.linkStyle = LinkStyle.Solid;
                addLinkParams.a = pA;
                addLinkParams.b = pB;
                addLinkCylinderMesh(addLinkOptions, addLinkParams);
                addLinkParams.a = pB;
                addLinkParams.b = pA;
                addLinkCylinderMesh(addLinkOptions, addLinkParams);
            }
        }
    }

    return Shape.create(
        'Interactions',
        interactions,
        MeshBuilder.getMesh(mesh),
        (g) => colors.get(g) ?? 0 as Color,
        (g) => 1,
        (g) => tooltips.get(g) ?? ''
    );
}


function cylinder(mesh: MeshBuilder.State, a: Vec3, b: Vec3, radius: number, style: 'dashed' | 'solid') {
    const props: BasicCylinderProps = {
        radiusBottom: radius,
        radiusTop: radius,
        topCap: true,
        bottomCap: true,
    };

    if (style === 'dashed') {
        const dist = Vec3.distance(a, b);
        const count = Math.ceil(dist / (2 * radius));
        addFixedCountDashedCylinder(mesh, a, b, 1.0, count, true, props);
    } else {
        if (style !== 'solid') {
            console.warn(`Unknown style '${style}', using 'solid' instead.`);
        }
        addSimpleCylinder(mesh, a, b, props);
    }
}