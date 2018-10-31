/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Color } from './color';

export interface BaseParam<T> {
    label: string
    description: string
    defaultValue: T
}

export interface ValueParam<T> extends BaseParam<T> {
    type: 'value'
}
export function ValueParam<T>(label: string, description: string, defaultValue: T): ValueParam<T> {
    return { type: 'value', label, description, defaultValue }
}

export interface SelectParam<T extends string> extends BaseParam<T> {
    type: 'select'
    /** array of (value, label) tupels */
    options: [T, string][]
}
export function SelectParam<T extends string>(label: string, description: string, defaultValue: T, options: [T, string][]): SelectParam<T> {
    return { type: 'select', label, description, defaultValue, options }
}

export interface MultiSelectParam<E extends string, T = E[]> extends BaseParam<T> {
    type: 'multi-select'
    /** array of (value, label) tupels */
    options: [E, string][]
}
export function MultiSelectParam<E extends string, T = E[]>(label: string, description: string, defaultValue: T, options: [E, string][]): MultiSelectParam<E, T> {
    return { type: 'multi-select', label, description, defaultValue, options }
}

export interface BooleanParam extends BaseParam<boolean> {
    type: 'boolean'
}
export function BooleanParam(label: string, description: string, defaultValue: boolean): BooleanParam {
    return { type: 'boolean', label, description, defaultValue }
}

export interface RangeParam extends BaseParam<number> {
    type: 'range'
    min: number
    max: number
    /** if an `integer` parse value with parseInt, otherwise use parseFloat */
    step: number
}
export function RangeParam(label: string, description: string, defaultValue: number, min: number, max: number, step: number): RangeParam {
    return { type: 'range', label, description, defaultValue, min, max, step }
}

export interface TextParam extends BaseParam<string> {
    type: 'text'
}
export function TextParam(label: string, description: string, defaultValue: string): TextParam {
    return { type: 'text', label, description, defaultValue }
}

export interface ColorParam extends BaseParam<Color> {
    type: 'color'
}
export function ColorParam(label: string, description: string, defaultValue: Color): ColorParam {
    return { type: 'color', label, description, defaultValue }
}

export interface NumberParam extends BaseParam<number> {
    type: 'number'
    min: number
    max: number
    /** if an `integer` parse value with parseInt, otherwise use parseFloat */
    step: number
}
export function NumberParam(label: string, description: string, defaultValue: number, min: number, max: number, step: number): NumberParam {
    return { type: 'number', label, description, defaultValue, min, max, step }
}

export type Param = ValueParam<any> | SelectParam<any> | MultiSelectParam<any> | BooleanParam | RangeParam | TextParam | ColorParam | NumberParam
export type Params = { [k: string]: Param }

export function paramDefaultValues<T extends Params>(params: T) {
    const d: { [k: string]: any } = {}
    Object.keys(params).forEach(k => d[k] = params[k].defaultValue)
    return d as { [k in keyof T]: T[k]['defaultValue'] }
}