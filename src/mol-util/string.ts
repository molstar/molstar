/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

const reLine = /^/mg
export function indentString(str: string, count: number, indent: string) {
    return count === 0 ? str : str.replace(reLine, indent.repeat(count))
}

/** Add space between camelCase text. */
export function splitCamelCase(str: string) {
    return str.replace(/([a-z\xE0-\xFF])([A-Z\xC0\xDF])/g, '$1 $2')
}

/** Split camelCase text and capitalize. */
export function camelCaseToWords(str: string) {
    return capitalize(splitCamelCase(str))
}

export const lowerCase = (str: string) => str.toLowerCase()
export const upperCase = (str: string) => str.toUpperCase()

/** Uppercase the first character of each word. */
export function capitalize(str: string) {
    return str.toLowerCase().replace(/^\w|\s\w/g, upperCase);
}

export function splitSnakeCase(str: string) {
    return str.replace(/_/g, ' ')
}

export function snakeCaseToWords(str: string) {
    return capitalize(splitSnakeCase(str))
}

export function stringToWords(str: string) {
    return capitalize(splitCamelCase(splitSnakeCase(str)))
}