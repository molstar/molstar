/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Color as ColorData } from './color';
import { shallowClone } from 'mol-util';
import { Vec2 } from 'mol-math/linear-algebra';
import { camelCaseToWords } from './string';

export namespace ParamDefinition {
    export interface Info {
        label?: string
        description?: string
    }

    export interface Base<T> extends Info {
        defaultValue: T
    }

    export interface Value<T> extends Base<T> {
        type: 'value'
    }
    export function Value<T>(defaultValue: T, info: Info = {}): Value<T> {
        return { type: 'value', defaultValue, ...info }
    }

    export interface Select<T extends string> extends Base<T> {
        type: 'select'
        /** array of (value, label) tuples */
        options: [T, string][]
    }
    export function Select<T extends string>(defaultValue: T, options: [T, string][], info: Info = {}): Select<T> {
        return { type: 'select', defaultValue, options, ...info }
    }

    export interface MultiSelect<E extends string, T = E[]> extends Base<T> {
        type: 'multi-select'
        /** array of (value, label) tuples */
        options: [E, string][]
    }
    export function MultiSelect<E extends string, T = E[]>(defaultValue: T, options: [E, string][], info: Info = {}): MultiSelect<E, T> {
        return { type: 'multi-select', defaultValue, options, ...info }
    }

    export interface Boolean extends Base<boolean> {
        type: 'boolean'
    }
    export function Boolean(defaultValue: boolean, info: Info = {}): Boolean {
        return { type: 'boolean', defaultValue, ...info }
    }

    export interface Text extends Base<string> {
        type: 'text'
    }
    export function Text(defaultValue: string = '', info: Info = {}): Text {
        return { type: 'text', defaultValue, ...info }
    }

    export interface Color extends Base<ColorData> {
        type: 'color'
    }
    export function Color(defaultValue: ColorData, info: Info = {}): Color {
        return { type: 'color', defaultValue, ...info }
    }

    export interface Numeric extends Base<number> {
        type: 'number'
        /** If given treat as a range. */
        min?: number
        /** If given treat as a range. */
        max?: number
        /**
         * If given treat as a range.
         * If an `integer` parse value with parseInt, otherwise use parseFloat.
         */
        step?: number
    }
    export function Numeric(defaultValue: number, range: { min?: number, max?: number, step?: number } = {}, info: Info = {}): Numeric {
        return { type: 'number', defaultValue, ...range, ...info }
    }

    export interface Interval extends Base<[number, number]> {
        type: 'interval'
    }
    export function Interval(defaultValue: [number, number], info: Info = {}): Interval {
        return { type: 'interval', defaultValue, ...info }
    }

    export interface LineGraph extends Base<Vec2[]> {
        type: 'line-graph'
    }
    export function LineGraph(defaultValue: Vec2[], info: Info = {}): LineGraph {
        return { type: 'line-graph', defaultValue, ...info }
    }

    export interface Group<T> extends Base<T> {
        type: 'group',
        params: Params
    }
    export function Group<P extends Params>(params: P, info: Info = {}): Group<Values<P>> {
        return { type: 'group', defaultValue: getDefaultValues(params) as any, params };
    }

    export interface Mapped<T> extends Base<{ name: string, params: T }> {
        type: 'mapped',
        select: Select<string>,
        map(name: string): Any
    }
    export function Mapped<T>(defaultKey: string, names: [string, string][], map: Mapped<T>['map'], info: Info = {}): Mapped<T> {
        return {
            type: 'mapped',
            defaultValue: { name: defaultKey, params: map(defaultKey).defaultValue as any },
            select: Select<string>(defaultKey, names, info),
            map
        };
    }

    export interface Converted<T, C> extends Base<T> {
        type: 'converted',
        convertedControl: Base<C>,
        fromValue(v: T): C,
        toValue(v: C): T
    }

    export type Any = Value<any> | Select<any> | MultiSelect<any> | Boolean | Text | Color | Numeric | Interval | LineGraph | Group<any> | Mapped<any> | Converted<any, any>

    export type Params = { [k: string]: Any }
    export type Values<T extends Params> = { [k in keyof T]: T[k]['defaultValue'] }

    export function getDefaultValues<T extends Params>(params: T) {
        const d: { [k: string]: any } = {}
        Object.keys(params).forEach(k => d[k] = params[k].defaultValue)
        return d as Values<T>
    }

    export function getLabels<T extends Params>(params: T) {
        const d: { [k: string]: string } = {}
        Object.keys(params).forEach(k => {
            const label = params[k].label
            d[k] = label === undefined ? camelCaseToWords(k) : label
        })
        return d as { [k in keyof T]: string }
    }

    export function clone<P extends Params>(params: P): P {
        return shallowClone(params)
    }

    /**
     * List of [error text, pathToValue]
     * i.e. ['Missing Nested Id', ['group1', 'id']]
     */
    export type ParamErrors = [string, string | string[]][]

    export function validate(params: Params, values: any): ParamErrors | undefined {
        // TODO
        return void 0;
    }
}