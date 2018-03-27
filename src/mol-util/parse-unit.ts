/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

/*
 * This code has been modified from https://github.com/mattdesl/parse-unit,
 * copyright (c) 2014 Matt DesLauriers. MIT License
 */

const reUnit = /[\d.\-\+]*\s*(.*)/

export default function parseUnit(str: string, out: [number, string] = [ 0, '' ]) {
    str = String(str)
    const num = parseFloat(str)
    out[0] = num
    const m = str.match(reUnit)
    if (m) out[1] = m[1] || ''
    return out
}