/**
 * Copyright (c) 2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Structure, StructureElement, StructureProperties, Unit } from 'mol-model/structure';
import { StateTransformer } from 'mol-state';
import { StructureLabels3D } from '../state/transforms/representation';
import { ShapeRepresentation } from 'mol-repr/shape/representation';
import { Vec3 } from 'mol-math/linear-algebra';
import { Text } from 'mol-geo/geometry/text/text';
import { TextBuilder } from 'mol-geo/geometry/text/text-builder';
import { Shape } from 'mol-model/shape';
import { ColorNames } from 'mol-util/color/tables';
import { RuntimeContext } from 'mol-task';
import { ParamDefinition as PD } from 'mol-util/param-definition';
import { BoundaryHelper } from 'mol-math/geometry/boundary-helper';

interface LabelsData {
    texts: string[],
    positions: Vec3[],
    sizes: number[],
    depths: number[]
}

function getLabelsText(data: LabelsData, props: PD.Values<Text.Params>, text?: Text) {
    const { texts, positions, depths } = data
    const textBuilder = TextBuilder.create(props, texts.length * 10, texts.length * 10 / 2, text)
    for (let i = 0, il = texts.length; i < il; ++i) {
        const p = positions[i]
        textBuilder.add(texts[i], p[0], p[1], p[2], depths[i], i)
    }
    return textBuilder.getText()
}

export async function getLabelRepresentation(ctx: RuntimeContext, structure: Structure, params: StateTransformer.Params<StructureLabels3D>, prev?: ShapeRepresentation<LabelsData, Text, Text.Params>) {
    const repr = prev || ShapeRepresentation(getLabelsShape, Text.Utils);
    const data = getLabelData(structure, params);
    await repr.createOrUpdate(params.options, data).runInContext(ctx);
    repr.setState({ pickable: false })
    return repr;
}

function getLabelsShape(ctx: RuntimeContext, data: LabelsData, props: PD.Values<Text.Params>, shape?: Shape<Text>) {
    const geo = getLabelsText(data, props, shape && shape.geometry);
    return Shape.create('Scene Labels', data, geo, () => ColorNames.dimgrey, g => data.sizes[g], () => '')
}

const boundaryHelper = new BoundaryHelper();
function getLabelData(structure: Structure, params: StateTransformer.Params<StructureLabels3D>): LabelsData {
    if (params.target.name === 'static-text') {
        return getLabelDataStatic(structure, params.target.params.value, params.target.params.size || 1, params.target.params.position || 'middle-center');
    } else {
        return getLabelDataComputed(structure, params.target.name);
    }

}

function getLabelDataStatic(structure: Structure, text: string, size: number, position: Text.Params['attachment']['defaultValue']): LabelsData {
    const boundary = structure.boundary.sphere;
    let oX = 0, oY = 0;
    if (position.indexOf('left') >= 0) oX = -boundary.radius;
    if (position.indexOf('right') >= 0) oX = boundary.radius;
    if (position.indexOf('top') >= 0) oY = boundary.radius;
    if (position.indexOf('bottom') >= 0) oY = -boundary.radius;
    return {
        texts: [text],
        positions: [Vec3.add(Vec3.zero(), boundary.center, Vec3.create(oX, oY, 0))],
        sizes: [size],
        depths: [boundary.radius + Math.sqrt(oX * oX + oY * oY)]
    };
}

function getLabelDataComputed(structure: Structure, level: 'elements' | 'residues'): LabelsData {
    const data: LabelsData = { texts: [], positions: [], sizes: [], depths: [] };

    const l = StructureElement.create();
    const { units } = structure;

    const { auth_atom_id } = StructureProperties.atom;
    const { auth_seq_id, auth_comp_id } = StructureProperties.residue;
    const { auth_asym_id } = StructureProperties.chain;
    const p = Vec3.zero();

    for (const unit of units) {
        // TODO: support coarse models

        if (unit.kind !== Unit.Kind.Atomic) continue;
        l.unit = unit;
        const elements = unit.elements;

        const pos = unit.conformation.position;

        if (level === 'elements') {
            for (let j = 0, _j = elements.length; j < _j; j++) {
                l.element = elements[j];

                pos(l.element, p);
                data.texts.push(auth_atom_id(l));
                data.positions.push(Vec3.clone(p));
                data.sizes.push(1);
                data.depths.push(2);
            }
        } else {
            const residueIndex = unit.model.atomicHierarchy.residueAtomSegments.index;

            let i = 0, len = elements.length;
            while (i < len) {
                const start = i, rI = residueIndex[elements[i]];
                i++;
                while (i < len && residueIndex[elements[i]] === rI) i++;

                boundaryHelper.reset(0);
                for (let eI = start; eI < i; eI++) {
                    pos(elements[eI], p);
                    boundaryHelper.boundaryStep(p, 0);
                }
                boundaryHelper.finishBoundaryStep();
                for (let eI = start; eI < i; eI++) {
                    pos(elements[eI], p);
                    boundaryHelper.extendStep(p, 0);
                }

                l.element = elements[start];

                data.texts.push(`${auth_comp_id(l)} ${auth_seq_id(l)}:${auth_asym_id(l)}`);
                data.positions.push(Vec3.clone(boundaryHelper.center));
                data.sizes.push(Math.max(1, boundaryHelper.radius / 5));
                data.depths.push(boundaryHelper.radius);
            }
        }
    }

    return data;
}