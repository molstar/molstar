/**
 * Copyright (c) 2017 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { OrderedSet, Iterator } from 'mol-data/int'
import { UniqueArray } from 'mol-data/generic'
import { SymmetryOperator } from 'mol-math/geometry/symmetry-operator'
import { Model, Format } from '../model'
import Unit from './unit'
import ElementSet from './element/set'
import ElementGroup from './element/group'
import Element from './element'

// A structure is a pair of "units" and an element set.
// Each unit contains the data and transformation of its corresponding elements.
interface Structure {
    readonly units: ReadonlyArray<Unit>,
    readonly elements: ElementSet
}

namespace Structure {
    export function create(units: ReadonlyArray<Unit>, elements: ElementSet): Structure { return { units, elements }; }
    export function Empty(units: ReadonlyArray<Unit>): Structure { return create(units, ElementSet.Empty); };

    export function ofData(format: Format) {
        const models = Model.create(format);
        return models.map(ofModel);
    }

    export function ofModel(model: Model): Structure {
        const chains = model.hierarchy.chainSegments;
        const builder = Builder();

        for (let c = 0; c < chains.count; c++) {
            const group = ElementGroup.createNew(OrderedSet.ofBounds(chains.segments[c], chains.segments[c + 1]));
            const unit = Unit.createAtomic(model, SymmetryOperator.Default, group);
            builder.add(unit, unit.fullGroup);
        }

        return builder.getStructure();
    }

    export interface Builder {
        add(unit: Unit, elements: ElementGroup): void,
        addUnit(unit: Unit): void,
        setElements(unitId: number, elements: ElementGroup): void,
        getStructure(): Structure,
        readonly elementCount: number
    }

    class BuilderImpl implements Builder {
        private _unitId = 0;
        private units: Unit[] = [];
        private elements = ElementSet.Generator();
        elementCount = 0;

        add(unit: Unit, elements: ElementGroup) { const id = this.addUnit(unit); this.setElements(id, elements); }
        addUnit(unit: Unit) { const id = this._unitId++; this.units[id] = unit; return id; }
        setElements(unitId: number, elements: ElementGroup) { this.elements.add(unitId, elements); this.elementCount += ElementGroup.size(elements); }
        getStructure(): Structure { return this.elementCount > 0 ? Structure.create(this.units, this.elements.getSet()) : Empty(this.units); }
    }

    export function Builder(): Builder { return new BuilderImpl(); }

    /** Transient = location gets overwritten when move() is called. */
    export function elementLocationsTransient(s: Structure): Iterator<Element.Location> {
        const l = Element.Location();
        const update = Element.updateLocation;
        return Iterator.map(ElementSet.elements(s.elements), a => update(s, l, a));
    }

    export function getModels(s: Structure) {
        const { units, elements } = s;
        const arr = UniqueArray.create<Model['id'], Model>();
        const ids = ElementSet.unitIds(elements);
        for (let i = 0; i < ids.length; i++) {
            const u = units[ids[i]];
            UniqueArray.add(arr, u.model.id, u.model);
        }
        return arr.array;
    }

    // TODO: "lift" atom set operators?
    // TODO: "diff"
}

export default Structure