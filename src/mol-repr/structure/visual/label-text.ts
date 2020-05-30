/**
 * Copyright (c) 2019-2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { ParamDefinition as PD } from '../../../mol-util/param-definition';
import { VisualUpdateState } from '../../../mol-repr/util';
import { VisualContext } from '../../../mol-repr/visual';
import { Structure, StructureElement, StructureProperties } from '../../../mol-model/structure';
import { Theme } from '../../../mol-theme/theme';
import { Text } from '../../../mol-geo/geometry/text/text';
import { TextBuilder } from '../../../mol-geo/geometry/text/text-builder';
import { ComplexTextVisual, ComplexTextParams, ComplexVisual } from '../complex-visual';
import { ElementIterator, getSerialElementLoci, eachSerialElement } from './util/element';
import { ColorNames } from '../../../mol-util/color/names';
import { Vec3 } from '../../../mol-math/linear-algebra';
import { PhysicalSizeTheme } from '../../../mol-theme/size/physical';
import { BoundaryHelper } from '../../../mol-math/geometry/boundary-helper';

export const LabelTextParams = {
    ...ComplexTextParams,
    background: PD.Boolean(true),
    backgroundMargin: PD.Numeric(0, { min: 0, max: 1, step: 0.01 }),
    backgroundColor: PD.Color(ColorNames.black),
    backgroundOpacity: PD.Numeric(0.5, { min: 0, max: 1, step: 0.01 }),
    level: PD.Select('residue', [['chain', 'Chain'], ['residue', 'Residue'], ['element', 'Element']] as const),
    chainScale: PD.Numeric(10, { min: 0, max: 20, step: 0.1 }),
    residueScale: PD.Numeric(1, { min: 0, max: 20, step: 0.1 }),
    elementScale: PD.Numeric(0.5, { min: 0, max: 20, step: 0.1 }),
};
export type LabelTextParams = typeof LabelTextParams
export type LabelTextProps = PD.Values<LabelTextParams>
export type LabelLevels = LabelTextProps['level']

export function LabelTextVisual(materialId: number): ComplexVisual<LabelTextParams> {
    return ComplexTextVisual<LabelTextParams>({
        defaultProps: PD.getDefaultValues(LabelTextParams),
        createGeometry: createLabelText,
        createLocationIterator: ElementIterator.fromStructure,
        getLoci: getSerialElementLoci,
        eachLocation: eachSerialElement,
        setUpdateState: (state: VisualUpdateState, newProps: PD.Values<LabelTextParams>, currentProps: PD.Values<LabelTextParams>) => {
            state.createGeometry = (
                newProps.level !== currentProps.level ||
                (newProps.level === 'chain' && newProps.chainScale !== currentProps.chainScale) ||
                (newProps.level === 'residue' && newProps.residueScale !== currentProps.residueScale) ||
                (newProps.level === 'element' && newProps.elementScale !== currentProps.elementScale)
            );
        }
    }, materialId);
}

function createLabelText(ctx: VisualContext, structure: Structure, theme: Theme, props: LabelTextProps, text?: Text): Text {

    switch (props.level) {
        case 'chain': return createChainText(ctx, structure, theme, props, text);
        case 'residue': return createResidueText(ctx, structure, theme, props, text);
        case 'element': return createElementText(ctx, structure, theme, props, text);
    }
}

//

const tmpVec = Vec3();
const boundaryHelper = new BoundaryHelper('98');

function createChainText(ctx: VisualContext, structure: Structure, theme: Theme, props: LabelTextProps, text?: Text): Text {
    const l = StructureElement.Location.create(structure);
    const { units, serialMapping } = structure;
    const { auth_asym_id, label_asym_id } = StructureProperties.chain;
    const { cumulativeUnitElementCount } = serialMapping;

    const count = units.length;
    const { chainScale } = props;
    const builder = TextBuilder.create(props, count, count / 2, text);

    for (let i = 0, il = units.length; i < il; ++i) {
        const unit = units[i];
        l.unit = unit;
        l.element = unit.elements[0];
        const { center, radius } = unit.lookup3d.boundary.sphere;
        Vec3.transformMat4(tmpVec, center, unit.conformation.operator.matrix);
        const authId = auth_asym_id(l);
        const labelId = label_asym_id(l);
        const text = authId === labelId ? labelId : `${labelId} [${authId}]`;
        builder.add(text, tmpVec[0], tmpVec[1], tmpVec[2], radius, chainScale, cumulativeUnitElementCount[i]);
    }

    return builder.getText();
}

function createResidueText(ctx: VisualContext, structure: Structure, theme: Theme, props: LabelTextProps, text?: Text): Text {
    const l = StructureElement.Location.create(structure);
    const { units, serialMapping } = structure;
    const { label_comp_id } = StructureProperties.atom;
    const { auth_seq_id } = StructureProperties.residue;
    const { cumulativeUnitElementCount } = serialMapping;

    const count = structure.polymerResidueCount * 2;
    const { residueScale } = props;
    const builder = TextBuilder.create(props, count, count / 2, text);

    for (let i = 0, il = units.length; i < il; ++i) {
        const unit = units[i];
        const pos = unit.conformation.position;
        const { elements } = unit;
        l.unit = unit;
        l.element = unit.elements[0];

        const residueIndex = unit.model.atomicHierarchy.residueAtomSegments.index;
        const groupOffset = cumulativeUnitElementCount[i];

        let j = 0, jl = elements.length;
        while (j < jl) {
            const start = j, rI = residueIndex[elements[j]];
            j++;
            while (j < jl && residueIndex[elements[j]] === rI) j++;

            boundaryHelper.reset();
            for (let eI = start; eI < j; eI++) {
                pos(elements[eI], tmpVec);
                boundaryHelper.includePosition(tmpVec);
            }
            boundaryHelper.finishedIncludeStep();
            for (let eI = start; eI < j; eI++) {
                pos(elements[eI], tmpVec);
                boundaryHelper.radiusPosition(tmpVec);
            }

            l.element = elements[start];

            const { center, radius } = boundaryHelper.getSphere();
            const authSeqId = auth_seq_id(l);
            const compId = label_comp_id(l);

            const text = `${compId} ${authSeqId}`;
            builder.add(text, center[0], center[1], center[2], radius, residueScale, groupOffset + start);
        }
    }

    return builder.getText();
}

function createElementText(ctx: VisualContext, structure: Structure, theme: Theme, props: LabelTextProps, text?: Text): Text {
    const l = StructureElement.Location.create(structure);
    const { units, serialMapping } = structure;
    const { label_atom_id, label_alt_id } = StructureProperties.atom;
    const { cumulativeUnitElementCount } = serialMapping;

    const sizeTheme = PhysicalSizeTheme({}, { scale: 1 });

    const count = structure.elementCount;
    const { elementScale } = props;
    const builder = TextBuilder.create(props, count, count / 2, text);

    for (let i = 0, il = units.length; i < il; ++i) {
        const unit = units[i];
        const pos = unit.conformation.position;
        const { elements } = unit;
        l.unit = unit;

        const groupOffset = cumulativeUnitElementCount[i];

        for (let j = 0, _j = elements.length; j < _j; j++) {
            l.element = elements[j];
            pos(l.element, tmpVec);
            const atomId = label_atom_id(l);
            const altId = label_alt_id(l);
            const text = altId ? `${atomId}%${altId}` : atomId;
            builder.add(text, tmpVec[0], tmpVec[1], tmpVec[2], sizeTheme.size(l), elementScale, groupOffset + j);
        }
    }

    return builder.getText();
}