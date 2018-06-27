/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Unit, StructureProperties, Element } from 'mol-model/structure';

import { StructureColorDataProps } from '.';
import { ColorData, createElementColor, createUniformColor } from '../../../util/color-data';
import { ColorScale } from 'mol-util/color';

function getAsymId(unit: Unit): Element.Property<string> {
    switch (unit.kind) {
        case Unit.Kind.Atomic:
            return StructureProperties.chain.label_asym_id
        case Unit.Kind.Spheres:
        case Unit.Kind.Gaussians:
            return StructureProperties.coarse.asym_id
    }
    throw new Error('unhandled unit kind')
}

export function chainIdColorData(props: StructureColorDataProps, locationFn: (l: Element.Location, renderElementIdx: number) => void, colorData?: ColorData) {
    const { group: { units }, elementCount } = props
    const unit = units[0]

    const map = unit.model.properties.asymIdSerialMap
    const count = map.size

    const domain = [ 0, count - 1 ]
    const scale = ColorScale.create({ domain })
    const asym_id = getAsymId(unit)

    const l = Element.Location()
    l.unit = unit

    return createElementColor({
        colorFn: (renderElementIdx: number) => {
            locationFn(l, renderElementIdx)
            return scale.color(map.get(asym_id(l)) || 0)
        },
        elementCount
    }, colorData)
}

export function chainIdElementColorData(props: StructureColorDataProps, colorData?: ColorData) {
    const elements = props.group.units[0].elements
    function locationFn(l: Element.Location, renderElementIdx: number) {
        l.element = elements[renderElementIdx]
    }
    return chainIdColorData(props, locationFn, colorData)
}

export function chainIdLinkColorData(props: StructureColorDataProps, colorData?: ColorData): ColorData {
    const unit = props.group.units[0]
    const elements = unit.elements
    let locationFn: (l: Element.Location, renderElementIdx: number) => void
    switch (unit.kind) {
        case Unit.Kind.Atomic:
            const { a } = unit.links
            locationFn = (l: Element.Location, renderElementIdx: number) => {
                l.element = elements[a[renderElementIdx]]
            }
            return chainIdColorData(props, locationFn, colorData)
        case Unit.Kind.Spheres:
        case Unit.Kind.Gaussians:
            // no chainId link color for coarse units, return uniform grey color
            return createUniformColor({ value: 0xCCCCCC }, colorData)
    }
}