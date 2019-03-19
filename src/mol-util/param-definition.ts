/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Color as ColorData } from './color';
import { shallowEqual } from 'mol-util';
import { Vec2 as Vec2Data, Vec3 as Vec3Data } from 'mol-math/linear-algebra';
import { deepClone } from './object';
import { Histogram as HistogramData } from 'mol-math/misc';

export namespace ParamDefinition {
    export interface Info {
        label?: string,
        description?: string,
        fieldLabels?: { [name: string]: string },
        isHidden?: boolean,
    }

    function setInfo<T extends Info>(param: T, info?: Info): T {
        if (!info) return param;
        if (info.label) param.label = info.label;
        if (info.description) param.description = info.description;
        if (info.fieldLabels) param.fieldLabels = info.fieldLabels;
        if (info.isHidden) param.isHidden = info.isHidden;
        return param;
    }

    export interface Base<T> extends Info {
        isOptional?: boolean,
        defaultValue: T
    }

    export interface Optional<T extends Any = Any> extends Base<T['defaultValue'] | undefined> {
        type: T['type']
    }

    export function Optional<T>(p: Base<T>): Base<T | undefined> {
        const ret = { ...p };
        ret.isOptional = true;
        return ret;
    }

    export interface Value<T> extends Base<T> {
        type: 'value'
    }
    export function Value<T>(defaultValue: T, info?: Info): Value<T> {
        return setInfo<Value<T>>({ type: 'value', defaultValue }, info);
    }

    export interface Select<T extends string | number> extends Base<T> {
        type: 'select'
        /** array of (value, label) tuples */
        options: [T, string][]
    }
    export function Select<T extends string | number>(defaultValue: T, options: [T, string][], info?: Info): Select<T> {
        return setInfo<Select<T>>({ type: 'select', defaultValue, options }, info)
    }

    export interface ColorScale<T extends string> extends Base<T> {
        type: 'color-scale'
        /** array of (value, label) tuples */
        options: [T, string][]
    }
    export function ColorScale<T extends string>(defaultValue: T, options: [T, string][], info?: Info): ColorScale<T> {
        return setInfo<ColorScale<T>>({ type: 'color-scale', defaultValue, options }, info)
    }

    export interface MultiSelect<E extends string, T = E[]> extends Base<T> {
        type: 'multi-select'
        /** array of (value, label) tuples */
        options: [E, string][]
    }
    export function MultiSelect<E extends string, T = E[]>(defaultValue: T, options: [E, string][], info?: Info): MultiSelect<E, T> {
        return setInfo<MultiSelect<E, T>>({ type: 'multi-select', defaultValue, options }, info)
    }

    export interface Boolean extends Base<boolean> {
        type: 'boolean'
    }
    export function Boolean(defaultValue: boolean, info?: Info): Boolean {
        return setInfo<Boolean>({ type: 'boolean', defaultValue }, info)
    }

    export interface Text<T extends string = string> extends Base<T> {
        type: 'text'
    }
    export function Text<T extends string = string>(defaultValue: string = '', info?: Info): Text<T> {
        return setInfo<Text<T>>({ type: 'text', defaultValue: defaultValue as any }, info)
    }

    export interface Color extends Base<ColorData> {
        type: 'color'
    }
    export function Color(defaultValue: ColorData, info?: Info): Color {
        return setInfo<Color>({ type: 'color', defaultValue }, info)
    }

    export interface Vec3 extends Base<Vec3Data> {
        type: 'vec3'
    }
    export function Vec3(defaultValue: Vec3Data, info?: Info): Vec3 {
        return setInfo<Vec3>({ type: 'vec3', defaultValue }, info)
    }

    export interface FileParam extends Base<File> {
        type: 'file',
        accept?: string
    }
    export function File(info?: Info & { accept?: string }): FileParam {
        const ret = setInfo<FileParam>({ type: 'file', defaultValue: void 0 as any }, info);
        if (info && info.accept) ret.accept = info.accept;
        return ret;
    }

    export interface Range {
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
    function setRange<T extends Numeric | Interval>(p: T, range?: { min?: number, max?: number, step?: number }) {
        if (!range) return p;
        if (typeof range.min !== 'undefined') p.min = range.min;
        if (typeof range.max !== 'undefined') p.max = range.max;
        if (typeof range.step !== 'undefined') p.step = range.step;
        return p;
    }

    export interface Numeric extends Base<number>, Range {
        type: 'number'
    }
    export function Numeric(defaultValue: number, range?: { min?: number, max?: number, step?: number }, info?: Info): Numeric {
        return setInfo<Numeric>(setRange({ type: 'number', defaultValue }, range), info)
    }

    export interface Interval extends Base<[number, number]>, Range {
        type: 'interval'
    }
    export function Interval(defaultValue: [number, number], range?: { min?: number, max?: number, step?: number }, info?: Info): Interval {
        return setInfo<Interval>(setRange({ type: 'interval', defaultValue }, range), info)
    }

    export interface LineGraph extends Base<Vec2Data[]> {
        type: 'line-graph'
    }
    export function LineGraph(defaultValue: Vec2Data[], info?: Info): LineGraph {
        return setInfo<LineGraph>({ type: 'line-graph', defaultValue }, info);
    }

    export interface Histogram extends Base<number> {
        type: 'histogram'
        histogram: HistogramData
    }
    
    export function Histogram(defaultValue: number, histogram: HistogramData, info?: Info): Histogram {
        return setInfo<Histogram>({ type: 'histogram', histogram, defaultValue}, info);
    }

    export interface Group<T> extends Base<T> {
        type: 'group',
        params: Params,
        isExpanded?: boolean,
        isFlat?: boolean
    }
    export function Group<T>(params: For<T>, info?: Info & { isExpanded?: boolean, isFlat?: boolean }): Group<Normalize<T>> {
        const ret = setInfo<Group<Normalize<T>>>({ type: 'group', defaultValue: getDefaultValues(params as any as Params) as any, params: params as any as Params }, info);
        if (info && info.isExpanded) ret.isExpanded = info.isExpanded;
        if (info && info.isFlat) ret.isFlat = info.isFlat;
        return ret;
    }

    export interface NamedParams<T = any, K = string> { name: K, params: T }
    export type NamedParamUnion<P extends Params, K = keyof P> = K extends any ? NamedParams<P[K]['defaultValue'], K> : never
    export interface Mapped<T extends NamedParams<any, any>> extends Base<T> {
        type: 'mapped',
        select: Select<string>,
        map(name: string): Any
    }
    export function Mapped<T>(defaultKey: string, names: [string, string][], map: (name: string) => Any, info?: Info): Mapped<NamedParams<T>> {
        return setInfo<Mapped<NamedParams<T>>>({
            type: 'mapped',
            defaultValue: { name: defaultKey, params: map(defaultKey).defaultValue as any },
            select: Select<string>(defaultKey, names, info),
            map
        }, info);
    }
    export function MappedStatic<C extends Params>(defaultKey: keyof C, map: C, info?: Info & { options?: [keyof C, string][] }): Mapped<NamedParamUnion<C>> {
        const options: [string, string][] = info && info.options
            ? info.options as [string, string][]
            : Object.keys(map).map(k => [k, k]) as [string, string][];
        return setInfo<Mapped<NamedParamUnion<C>>>({
            type: 'mapped',
            defaultValue: { name: defaultKey, params: map[defaultKey].defaultValue } as any,
            select: Select<string>(defaultKey as string, options, info),
            map: key => map[key]
        }, info);
    }

    export interface ObjectList<T = any> extends Base<T[]> {
        type: 'object-list',
        element: Params,
        ctor(): T,
        getLabel(t: T): string
    }
    export function ObjectList<T>(element: For<T>, getLabel: (e: T) => string, info?: Info & { defaultValue?: T[], ctor?: () => T }): ObjectList<Normalize<T>> {
        return setInfo<ObjectList<Normalize<T>>>({ type: 'object-list', element: element as any as Params, getLabel, ctor: _defaultObjectListCtor, defaultValue: (info && info.defaultValue) || []  });
    }
    function _defaultObjectListCtor(this: ObjectList) { return getDefaultValues(this.element) as any; }

    export interface Converted<T, C> extends Base<T> {
        type: 'converted',
        converted: Any,
        /** converts from prop value to display value */
        fromValue(v: T): C,
        /** converts from display value to prop value */
        toValue(v: C): T
    }
    export function Converted<T, C extends Any>(fromValue: (v: T) => C['defaultValue'], toValue: (v: C['defaultValue']) => T, converted: C): Converted<T, C['defaultValue']> {
        return { type: 'converted', defaultValue: toValue(converted.defaultValue), converted, fromValue, toValue };
    }

    export interface Conditioned<T, P extends Base<T>, C = { [k: string]: P }> extends Base<T> {
        type: 'conditioned',
        select: Select<string>,
        conditionParams: C
        conditionForValue(v: T): keyof C
        conditionedValue(v: T, condition: keyof C): T,
    }
    export function Conditioned<T, P extends Base<T>, C = { [k: string]: P }>(defaultValue: T, conditionParams: C, conditionForValue: (v: T) => keyof C, conditionedValue: (v: T, condition: keyof C) => T): Conditioned<T, P, C> {
        const options = Object.keys(conditionParams).map(k => [k, k]) as [string, string][];
        return { type: 'conditioned', select: Select<string>(conditionForValue(defaultValue) as string, options), defaultValue, conditionParams, conditionForValue, conditionedValue };
    }

    export interface ScriptExpression extends Base<{ language: 'mol-script', expression: string }> {
        type: 'script-expression'
    }
    export function ScriptExpression(defaultValue: ScriptExpression['defaultValue'], info?: Info): ScriptExpression {
        return setInfo<ScriptExpression>({ type: 'script-expression', defaultValue }, info)
    }

    export type Any =
        | Value<any> | Select<any> | MultiSelect<any> | Boolean | Text | Color | Vec3 | Numeric | FileParam | Interval | LineGraph | Histogram
        | ColorScale<any> | Group<any> | Mapped<any> | Converted<any, any> | Conditioned<any, any, any> | ScriptExpression | ObjectList

    export type Params = { [k: string]: Any }
    export type Values<T extends Params> = { [k in keyof T]: T[k]['defaultValue'] }

    type Optionals<P> = { [K in keyof P]-?: undefined extends P[K] ? K : never }[keyof P]
    type NonOptionals<P> = { [K in keyof P]-?: undefined extends P[K] ? never: K }[keyof P]
    export type Normalize<P> = Pick<P, NonOptionals<P>> & Partial<Pick<P, Optionals<P>>>
    export type For<P> = { [K in keyof P]-?: Base<P[K]> }

    export function getDefaultValues<T extends Params>(params: T) {
        const d: { [k: string]: any } = {}
        for (const k of Object.keys(params)) {
            if (params[k].isOptional) continue;
            d[k] = params[k].defaultValue;
        }
        return d as Values<T>;
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
        } else if (p.type === 'multi-select') {
            const u = a as MultiSelect<any>['defaultValue'], v = b as MultiSelect<any>['defaultValue'];
            if (u.length !== v.length) return false;
            if (u.length < 10) {
                for (let i = 0, _i = u.length; i < _i; i++) {
                    if (u[i] === v[i]) continue;
                    if (v.indexOf(u[i]) < 0) return false;
                }
            } else {
                // TODO: should the value of multiselect be a set?
                const vSet = new Set(v);
                for (let i = 0, _i = u.length; i < _i; i++) {
                    if (u[i] === v[i]) continue;
                    if (!vSet.has(u[i])) return false;
                }
            }
            return true;
        } else if (p.type === 'interval') {
            return a[0] === b[0] && a[1] === b[1];
        } else if (p.type === 'line-graph') {
            const u = a as LineGraph['defaultValue'], v = b as LineGraph['defaultValue'];
            if (u.length !== v.length) return false;
            for (let i = 0, _i = u.length; i < _i; i++) {
                if (!Vec2Data.areEqual(u[i], v[i])) return false;
            }
            return true;
        } else if (p.type === 'vec3') {
            return Vec3Data.equals(a, b);
        } else if (p.type === 'script-expression') {
            const u = a as ScriptExpression['defaultValue'], v = b as ScriptExpression['defaultValue'];
            return u.language === v.language && u.expression === v.expression;
        } else if (p.type === 'object-list') {
            const u = a as ObjectList['defaultValue'], v = b as ObjectList['defaultValue'];
            const l = u.length;
            if (l !== v.length) return false;
            for (let i = 0; i < l; i++) {
                if (!areEqual(p.element, u[i], v[i])) return false;
            }
            return true;
        } else if (typeof a === 'object' && typeof b === 'object') {
            return shallowEqual(a, b);
        }

        // a === b was checked at the top.
        return false;
    }
}