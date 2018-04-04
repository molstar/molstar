/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

/*
 * This code has been modified from https://github.com/mikolalysenko/to-px,
 * copyright (c) 2015 Mikola Lysenko. MIT License
 */

import parseUnit from './parse-unit'

const PIXELS_PER_INCH = 96

function getPropertyInPX(element: Element, prop: string) {
    const parts = parseUnit(getComputedStyle(element).getPropertyValue(prop))
    return parts[0] * toPixels(parts[1], element)
}

// This brutal hack is needed
function getSizeBrutal(unit: string, element: Element) {
    const testDIV = document.createElement('div')
    testDIV.style.setProperty('font-size', '128' + unit)
    element.appendChild(testDIV)
    const size = getPropertyInPX(testDIV, 'font-size') / 128
    element.removeChild(testDIV)
    return size
}

export default function toPixels(str: string, element: Element = document.body): number {
    str = (str || 'px').trim().toLowerCase()
    switch (str) {
        case '%':  // Ambiguous, not sure if we should use width or height
            return element.clientHeight / 100.0
        case 'ch':
        case 'ex':
            return getSizeBrutal(str, element)
        case 'em':
            return getPropertyInPX(element, 'font-size')
        case 'rem':
            return getPropertyInPX(document.body, 'font-size')
        case 'vw':
            return window.innerWidth/100
        case 'vh':
            return window.innerHeight/100
        case 'vmin':
            return Math.min(window.innerWidth, window.innerHeight) / 100
        case 'vmax':
            return Math.max(window.innerWidth, window.innerHeight) / 100
        case 'in':
            return PIXELS_PER_INCH
        case 'cm':
            return PIXELS_PER_INCH / 2.54
        case 'mm':
            return PIXELS_PER_INCH / 25.4
        case 'pt':
            return PIXELS_PER_INCH / 72
        case 'pc':
            return PIXELS_PER_INCH / 6
  }
  return 1
}