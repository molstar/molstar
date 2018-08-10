/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Iterator } from 'mol-data';
import { Unit, StructureElement, Structure, Link } from 'mol-model/structure';
import { NullLocation, Location } from 'mol-model/location';

export interface LocationValue {
    location: Location
    index: number
    elementIndex: number
    instanceIndex: number
    isSecondary: boolean
}

export const NullLocationValue: LocationValue = {
    location: NullLocation,
    index: 0,
    elementIndex: 0,
    instanceIndex: 0,
    isSecondary: false
}

export interface LocationIterator extends Iterator<LocationValue> {
    readonly hasNext: boolean
    readonly isNextNewInstance: boolean
    readonly elementCount: number
    readonly instanceCount: number
    move(): LocationValue
    skipInstance(): void
}

type LocationGetter = (elementIndex: number, instanceIndex: number) => Location
type IsSecondaryGetter = (elementIndex: number, instanceIndex: number) => boolean

export function LocationIterator(elementCount: number, instanceCount: number, getLocation: LocationGetter, isSecondary: IsSecondaryGetter = () => false): LocationIterator {
    const value: LocationValue = {
        location: NullLocation as Location,
        index: 0,
        elementIndex: 0,
        instanceIndex: 0,
        isSecondary: false
    }

    let hasNext = value.elementIndex < elementCount
    let isNextNewInstance = false
    let elementIndex = 0
    let instanceIndex = 0

    return {
        get hasNext () { return hasNext },
        get isNextNewInstance () { return isNextNewInstance },
        get elementCount () { return elementCount },
        get instanceCount () { return instanceCount },
        move() {
            if (hasNext) {
                value.elementIndex = elementIndex
                value.instanceIndex = instanceIndex
                value.index = instanceIndex * elementCount + elementIndex
                value.location = getLocation(elementIndex, instanceIndex)
                value.isSecondary = isSecondary(elementIndex, instanceIndex)
                ++elementIndex
                if (elementIndex === elementCount) {
                    ++instanceIndex
                    isNextNewInstance = true
                    if (instanceIndex < instanceCount) elementIndex = 0
                } else {
                    isNextNewInstance = false
                }
                hasNext = elementIndex < elementCount
            }
            return value
        },
        skipInstance() {
            if (hasNext && value.instanceIndex === instanceIndex) {
                ++instanceIndex
                elementIndex = 0
                hasNext = instanceIndex < instanceCount
            }
        }
    }
}

export namespace StructureElementIterator {
    export function fromGroup(group: Unit.SymmetryGroup): LocationIterator {
        const unit = group.units[0]
        const elementCount = group.elements.length
        const instanceCount = group.units.length
        const location = StructureElement.create(unit)
        const getLocation = (elementIndex: number, instanceIndex: number) => {
            location.element = unit.elements[elementIndex]
            return location
        }
        return LocationIterator(elementCount, instanceCount, getLocation)
    }
}

export namespace LinkIterator {
    export function fromGroup(group: Unit.SymmetryGroup): LocationIterator {
        const unit = group.units[0]
        const elementCount = Unit.isAtomic(unit) ? unit.links.edgeCount * 2 : 0
        const instanceCount = group.units.length
        const location = StructureElement.create(unit)
        const getLocation = (elementIndex: number, instanceIndex: number) => {
            location.element = unit.elements[(unit as Unit.Atomic).links.a[elementIndex]]
            return location
        }
        return LocationIterator(elementCount, instanceCount, getLocation)
    }

    export function fromStructure(structure: Structure): LocationIterator {
        const elementCount = structure.links.bondCount
        const instanceCount = 1
        const location = Link.Location()
        const getLocation = (elementIndex: number, instanceIndex: number) => {
            const bond = structure.links.bonds[elementIndex]
            location.aUnit = bond.unitA
            location.aIndex = bond.indexA as StructureElement.UnitIndex
            location.bUnit = bond.unitB
            location.bIndex = bond.indexB as StructureElement.UnitIndex
            return location
        }
        return LocationIterator(elementCount, instanceCount, getLocation)
    }
}