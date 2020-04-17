/**
 * Copyright (c) 2019-2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { CustomElementProperty } from '../../mol-model-props/common/custom-element-property';
import { Model, ElementIndex } from '../../mol-model/structure';
import { Color } from '../../mol-util/color';

export const StripedResidues = CustomElementProperty.create<number>({
    label: 'Residue Stripes',
    name: 'basic-wrapper-residue-striping',
    getData(model: Model) {
        const map = new Map<ElementIndex, number>();
        const residueIndex = model.atomicHierarchy.residueAtomSegments.index;
        for (let i = 0, _i = model.atomicHierarchy.atoms._rowCount; i < _i; i++) {
            map.set(i as ElementIndex, residueIndex[i] % 2);
        }
        return { value: map };
    },
    coloring: {
        getColor(e) { return e === 0 ? Color(0xff0000) : Color(0x0000ff); },
        defaultColor: Color(0x777777)
    },
    getLabel(e) {
        return e === 0 ? 'Odd stripe' : 'Even stripe';
    }
});