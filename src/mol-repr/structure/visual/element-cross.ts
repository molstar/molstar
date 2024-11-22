/**
 * Copyright (c) 2021-2024 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { ParamDefinition as PD } from '../../../mol-util/param-definition';
import { UnitsVisual, UnitsLinesParams, UnitsLinesVisual } from '../units-visual';
import { VisualContext } from '../../visual';
import { Unit, Structure, StructureElement } from '../../../mol-model/structure';
import { Theme } from '../../../mol-theme/theme';
import { Vec3 } from '../../../mol-math/linear-algebra';
import { ElementIterator, getElementLoci, eachElement, makeElementIgnoreTest, getSerialElementLoci, eachSerialElement } from './util/element';
import { VisualUpdateState } from '../../util';
import { Sphere3D } from '../../../mol-math/geometry';
import { Lines } from '../../../mol-geo/geometry/lines/lines';
import { LinesBuilder } from '../../../mol-geo/geometry/lines/lines-builder';
import { bondCount } from '../../../mol-model-props/computed/chemistry/util';
import { hasUnitVisibleBonds } from './util/bond';
import { ComplexLinesParams, ComplexLinesVisual, ComplexVisual } from '../complex-visual';

// avoiding namespace lookup improved performance in Chrome (Aug 2020)
const v3add = Vec3.add;
const v3scaleAndAdd = Vec3.scaleAndAdd;
const v3unitX = Vec3.unitX;
const v3unitY = Vec3.unitY;
const v3unitZ = Vec3.unitZ;

export const ElementCrossParams = {
    ...UnitsLinesParams,
    lineSizeAttenuation: PD.Boolean(false),
    ignoreHydrogens: PD.Boolean(false),
    ignoreHydrogensVariant: PD.Select('all', PD.arrayToOptions(['all', 'non-polar'] as const)),
    traceOnly: PD.Boolean(false),
    crosses: PD.Select('lone', PD.arrayToOptions(['lone', 'all'] as const)),
    crossSize: PD.Numeric(0.35, { min: 0, max: 2, step: 0.01 }),
};
export type ElementCrossParams = typeof ElementCrossParams

export function createElementCross(ctx: VisualContext, unit: Unit, structure: Structure, theme: Theme, props: PD.Values<ElementCrossParams>, lines: Lines) {
    const { child } = structure;
    if (child && !child.unitMap.get(unit.id)) return Lines.createEmpty(lines);

    const elements = unit.elements;
    const n = elements.length;
    const builder = LinesBuilder.create(n, n / 10, lines);

    const p = Vec3();
    const s = Vec3();
    const e = Vec3();

    const c = unit.conformation;
    const ignore = makeElementIgnoreTest(structure, unit, props);

    const r = props.crossSize / 2;
    const lone = props.crosses === 'lone';

    const center = Vec3();
    let count = 0;

    for (let i = 0 as StructureElement.UnitIndex; i < n; ++i) {
        if (ignore && ignore(elements[i])) continue;
        if (lone && Unit.isAtomic(unit) && hasUnitVisibleBonds(unit, props) && bondCount(structure, unit, i) !== 0) continue;

        c.invariantPosition(elements[i], p);
        v3add(center, center, p);
        count += 1;

        v3scaleAndAdd(s, p, v3unitX, r);
        v3scaleAndAdd(e, p, v3unitX, -r);
        builder.add(s[0], s[1], s[2], e[0], e[1], e[2], i);
        v3scaleAndAdd(s, p, v3unitY, r);
        v3scaleAndAdd(e, p, v3unitY, -r);
        builder.add(s[0], s[1], s[2], e[0], e[1], e[2], i);
        v3scaleAndAdd(s, p, v3unitZ, r);
        v3scaleAndAdd(e, p, v3unitZ, -r);
        builder.add(s[0], s[1], s[2], e[0], e[1], e[2], i);
    }

    const l = builder.getLines();
    if (count === 0) return l;

    // re-use boundingSphere if it has not changed much
    let boundingSphere: Sphere3D;
    Vec3.scale(center, center, 1 / count);
    const oldBoundingSphere = lines ? Sphere3D.clone(lines.boundingSphere) : undefined;
    if (oldBoundingSphere && Vec3.distance(center, oldBoundingSphere.center) / oldBoundingSphere.radius < 0.1) {
        boundingSphere = oldBoundingSphere;
    } else {
        boundingSphere = Sphere3D.expand(Sphere3D(), unit.boundary.sphere, 1 * props.sizeFactor);
    }
    l.setBoundingSphere(boundingSphere);

    return l;
}

export function ElementCrossVisual(materialId: number): UnitsVisual<ElementCrossParams> {
    return UnitsLinesVisual<ElementCrossParams>({
        defaultProps: PD.getDefaultValues(ElementCrossParams),
        createGeometry: createElementCross,
        createLocationIterator: ElementIterator.fromGroup,
        getLoci: getElementLoci,
        eachLocation: eachElement,
        setUpdateState: (state: VisualUpdateState, newProps: PD.Values<ElementCrossParams>, currentProps: PD.Values<ElementCrossParams>) => {
            state.createGeometry = (
                newProps.ignoreHydrogens !== currentProps.ignoreHydrogens ||
                newProps.ignoreHydrogensVariant !== currentProps.ignoreHydrogensVariant ||
                newProps.traceOnly !== currentProps.traceOnly ||
                newProps.crosses !== currentProps.crosses ||
                newProps.crossSize !== currentProps.crossSize
            );
        }
    }, materialId);
}

//

export function createStructureElementCross(ctx: VisualContext, structure: Structure, theme: Theme, props: PD.Values<StructureElementCrossParams>, lines?: Lines): Lines {
    const { child } = structure;

    const { getSerialIndex } = structure.serialMapping;
    const structureElementCount = structure.elementCount;
    const builder = LinesBuilder.create(structureElementCount, structureElementCount / 2, lines);

    const p = Vec3();
    const s = Vec3();
    const e = Vec3();

    const r = props.crossSize / 2;
    const lone = props.crosses === 'lone';

    const center = Vec3();
    let count = 0;

    for (const unit of structure.units) {
        const childUnit = child?.unitMap.get(unit.id);
        if (child && !childUnit) return Lines.createEmpty(lines);

        const { elements, conformation: c } = unit;
        const elementCount = elements.length;

        const ignore = makeElementIgnoreTest(structure, unit, props);

        for (let i = 0 as StructureElement.UnitIndex; i < elementCount; i++) {
            if (ignore && ignore(elements[i])) continue;
            if (lone && Unit.isAtomic(unit) && hasUnitVisibleBonds(unit, props) && bondCount(structure, unit, i) !== 0) continue;

            c.position(elements[i], p);
            v3add(center, center, p);
            count += 1;

            const si = getSerialIndex(unit, elements[i]);
            v3scaleAndAdd(s, p, v3unitX, r);
            v3scaleAndAdd(e, p, v3unitX, -r);
            builder.add(s[0], s[1], s[2], e[0], e[1], e[2], si);
            v3scaleAndAdd(s, p, v3unitY, r);
            v3scaleAndAdd(e, p, v3unitY, -r);
            builder.add(s[0], s[1], s[2], e[0], e[1], e[2], si);
            v3scaleAndAdd(s, p, v3unitZ, r);
            v3scaleAndAdd(e, p, v3unitZ, -r);
            builder.add(s[0], s[1], s[2], e[0], e[1], e[2], si);
        }
    }

    const l = builder.getLines();
    if (count === 0) return l;

    // re-use boundingSphere if it has not changed much
    let boundingSphere: Sphere3D;
    Vec3.scale(center, center, 1 / count);
    const oldBoundingSphere = lines ? Sphere3D.clone(lines.boundingSphere) : undefined;
    if (oldBoundingSphere && Vec3.distance(center, oldBoundingSphere.center) / oldBoundingSphere.radius < 1.0) {
        boundingSphere = oldBoundingSphere;
    } else {
        boundingSphere = Sphere3D.expand(Sphere3D(), (child ?? structure).boundary.sphere, 1 * props.sizeFactor);
    }
    l.setBoundingSphere(boundingSphere);

    return l;
}

export const StructureElementCrossParams = {
    ...ComplexLinesParams,
    lineSizeAttenuation: PD.Boolean(false),
    ignoreHydrogens: PD.Boolean(false),
    ignoreHydrogensVariant: PD.Select('all', PD.arrayToOptions(['all', 'non-polar'] as const)),
    traceOnly: PD.Boolean(false),
    crosses: PD.Select('lone', PD.arrayToOptions(['lone', 'all'] as const)),
    crossSize: PD.Numeric(0.35, { min: 0, max: 2, step: 0.01 }),
};
export type StructureElementCrossParams = typeof StructureElementCrossParams

export function StructureElementCrossVisual(materialId: number): ComplexVisual<StructureElementCrossParams> {
    return ComplexLinesVisual<StructureElementCrossParams>({
        defaultProps: PD.getDefaultValues(StructureElementCrossParams),
        createGeometry: createStructureElementCross,
        createLocationIterator: ElementIterator.fromStructure,
        getLoci: getSerialElementLoci,
        eachLocation: eachSerialElement,
        setUpdateState: (state: VisualUpdateState, newProps: PD.Values<StructureElementCrossParams>, currentProps: PD.Values<StructureElementCrossParams>) => {
            state.createGeometry = (
                newProps.ignoreHydrogens !== currentProps.ignoreHydrogens ||
                newProps.ignoreHydrogensVariant !== currentProps.ignoreHydrogensVariant ||
                newProps.traceOnly !== currentProps.traceOnly ||
                newProps.crosses !== currentProps.crosses ||
                newProps.crossSize !== currentProps.crossSize
            );
        }
    }, materialId);
}
