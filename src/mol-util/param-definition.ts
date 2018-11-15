/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Color as ColorData } from './color';
import { shallowClone, shallowEqual } from 'mol-util';
import { Vec2 } from 'mol-math/linear-algebra';
import { deepClone } from './object';

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
        return { type: 'group', defaultValue: getDefaultValues(params) as any, params, ...info };
    }

    export interface NamedParams<T = any> { name: string, params: T }
    export interface Mapped<T> extends Base<NamedParams<T>> {
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

    export interface Converted<T, C> extends Base<C> {
        type: 'converted',
        param: Base<C>,
        /** converts from prop value to display value */
        fromValue(v: T): C,
        /** converts from display value to prop value */
        toValue(v: C): T
    }
    export function Converted<T, C>(param: Base<C>, fromValue: Converted<T, C>['fromValue'], toValue: Converted<T, C>['toValue']): Converted<T, C> {
        return { type: 'converted', param, defaultValue: param.defaultValue, fromValue, toValue, label: param.label, description: param.description };
    }

    export type Any = Value<any> | Select<any> | MultiSelect<any> | Boolean | Text | Color | Numeric | Interval | LineGraph | Group<any> | Mapped<any> | Converted<any, any>

    export type Params = { [k: string]: Any }
    export type Values<T extends Params> = { [k in keyof T]: T[k]['defaultValue'] }

    export function getDefaultValues<T extends Params>(params: T) {
        const d: { [k: string]: any } = {}
        Object.keys(params).forEach(k => d[k] = params[k].defaultValue)
        return d as Values<T>
    }

    export function clone<P extends Params>(params: P): P {
        return deepClone(params)
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

    export function areEqual(params: Params, a: any, b: any): boolean {
        if (a === b) return true;
        if (!a) return !b;
        if (!b) return !a;

        if (typeof a !== 'object' || typeof b !== 'object') return false;
        for (const k of Object.keys(params)) {
            if (!isParamEqual(params[k], a[k], b[k])) return false;
        }
        return true;
    }

    function isParamEqual(p: Any, a: any, b: any): boolean {
        if (a === b) return true;
        if (!a) return !b;
        if (!b) return !a;

        if (p.type === 'group') {
            return areEqual(p.params, a, b);
        } else if (p.type === 'mapped') {
            const u = a as NamedParams, v = b as NamedParams;
            if (!u) return !v;
            if (!u || !v) return false;
            if (u.name !== v.name) return false;
            const map = p.map(u.name);
            return isParamEqual(map, u.params, v.params);
        } else if (p.type === 'interval') {
            return a[0] === b[0] && a[1] === b[1];
        } else if (p.type === 'line-graph') {
            const u = a as LineGraph['defaultValue'], v = b as LineGraph['defaultValue'];
            if (u.length !== v.length) return false;
            for (let i = 0, _i = u.length; i < _i; i++) {
                if (!Vec2.areEqual(u[i], v[i])) return false;
            }
            return true;
        } else if (typeof a === 'object' && typeof b === 'object') {
            return shallowEqual(a, b);
        }

        // a === b was checked at the top.
        return false;
    }
}