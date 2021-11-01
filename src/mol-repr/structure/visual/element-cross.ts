/**
 * Copyright (c) 2021 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { ParamDefinition as PD } from '../../../mol-util/param-definition';
import { UnitsVisual, UnitsLinesParams, UnitsLinesVisual } from '../units-visual';
import { VisualContext } from '../../visual';
import { Unit, Structure, StructureElement } from '../../../mol-model/structure';
import { Theme } from '../../../mol-theme/theme';
import { Vec3 } from '../../../mol-math/linear-algebra';
import { ElementIterator, getElementLoci, eachElement, makeElementIgnoreTest } from './util/element';
import { VisualUpdateState } from '../../util';
import { Sphere3D } from '../../../mol-math/geometry';
import { Lines } from '../../../mol-geo/geometry/lines/lines';
import { LinesBuilder } from '../../../mol-geo/geometry/lines/lines-builder';
import { bondCount } from '../../../mol-model-props/computed/chemistry/util';

// avoiding namespace lookup improved performance in Chrome (Aug 2020)
const v3scaleAndAdd = Vec3.scaleAndAdd;
const v3unitX = Vec3.unitX;
const v3unitY = Vec3.unitY;
const v3unitZ = Vec3.unitZ;

export const ElementCrossParams = {
    ...UnitsLinesParams,
    lineSizeAttenuation: PD.Boolean(false),
    ignoreHydrogens: PD.Boolean(false),
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

    const pos = unit.conformation.invariantPosition;
    const ignore = makeElementIgnoreTest(structure, unit, props);

    const r = props.crossSize / 2;
    const lone = props.crosses === 'lone';

    for (let i = 0 as StructureElement.UnitIndex; i < n; ++i) {
        if (ignore && ignore(elements[i])) continue;
        if (lone && Unit.isAtomic(unit) && bondCount(structure, unit, i) !== 0) continue;

        pos(elements[i], p);
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

    const sphere = Sphere3D.expand(Sphere3D(), unit.boundary.sphere, 1 * props.sizeFactor);
    l.setBoundingSphere(sphere);

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
                newProps.traceOnly !== currentProps.traceOnly ||
                newProps.crosses !== currentProps.crosses ||
                newProps.crossSize !== currentProps.crossSize
            );
        }
    }, materialId);
}