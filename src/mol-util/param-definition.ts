/**
 * Copyright (c) 2018-2024 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Color as ColorData } from './color';
import { shallowEqualObjects } from './index';
import { Vec2 as Vec2Data, Vec3 as Vec3Data, Mat4 as Mat4Data, EPSILON } from '../mol-math/linear-algebra';
import { deepClone } from './object';
import { Script as ScriptData } from '../mol-script/script';
import { Legend } from './legend';
import { stringToWords } from './string';
import { getColorListFromName, ColorListName } from './color/lists';
import { Asset } from './assets';
import { ColorListEntry } from './color/color';

export namespace ParamDefinition {
    export interface Info {
        label?: string,
        description?: string,
        legend?: Legend,
        fieldLabels?: { [name: string]: string },
        isHidden?: boolean,
        shortLabel?: boolean,
        twoColumns?: boolean,
        isEssential?: boolean,
        category?: string,

        hideIf?: (currentGroup: any) => boolean,
        help?: (value: any) => { description?: string, legend?: Legend }
    }

    export const Essential = { isEssential: true };

    function setInfo<T extends Base<any>>(param: T, info?: Info): T {
        if (!info) return param;
        if (info.label) param.label = info.label;
        if (info.description) param.description = info.description;
        if (info.legend) param.legend = info.legend;
        if (info.fieldLabels) param.fieldLabels = info.fieldLabels;
        if (info.isHidden) param.isHidden = info.isHidden;
        if (info.shortLabel) param.shortLabel = info.shortLabel;
        if (info.twoColumns) param.twoColumns = info.twoColumns;
        if (info.isEssential) param.isEssential = info.isEssential;
        if (info.category) param.category = info.category;

        if (info.hideIf) param.hideIf = info.hideIf;
        if (info.help) param.help = info.help;
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

    export interface Select<T> extends Base<T> {
        type: 'select'
        /** array of (value, label) tuples */
        options: readonly (readonly [T, string] | readonly [T, string, string | undefined])[]
        cycle?: boolean
    }
    export function Select<T>(defaultValue: T, options: readonly (readonly [T, string] | readonly [T, string, string | undefined])[], info?: Info & { cycle?: boolean }): Select<T> {
        return setInfo<Select<T>>({ type: 'select', defaultValue: checkDefaultKey(defaultValue, options), options, cycle: info?.cycle }, info);
    }

    export interface MultiSelect<E extends string> extends Base<E[]> {
        type: 'multi-select'
        /** array of (value, label) tuples */
        options: readonly (readonly [E, string])[],
        emptyValue?: string
    }
    export function MultiSelect<E extends string>(defaultValue: E[], options: readonly (readonly [E, string])[], info?: Info & { emptyValue?: string }): MultiSelect<E> {
        // TODO: check if default value is a subset of options?
        const ret = setInfo<MultiSelect<E>>({ type: 'multi-select', defaultValue, options }, info);
        if (info?.emptyValue) ret.emptyValue = info.emptyValue;
        return ret;
    }

    export interface BooleanParam extends Base<boolean> {
        type: 'boolean'
    }
    export function Boolean(defaultValue: boolean, info?: Info): BooleanParam {
        return setInfo<BooleanParam>({ type: 'boolean', defaultValue }, info);
    }

    export interface Text<T extends string = string> extends Base<T> {
        type: 'text',
        multiline?: boolean,
        placeholder?: string,
        disableInteractiveUpdates?: boolean
    }
    export function Text<T extends string = string>(defaultValue: string = '', info?: Info & { multiline?: boolean, placeholder?: string, disableInteractiveUpdates?: boolean }): Text<T> {
        return setInfo<Text<T>>({ type: 'text', defaultValue: defaultValue as any, multiline: info?.multiline, placeholder: info?.placeholder, disableInteractiveUpdates: info?.disableInteractiveUpdates }, info);
    }

    export interface Color extends Base<ColorData> {
        type: 'color'
        isExpanded?: boolean
    }
    export function Color(defaultValue: ColorData, info?: Info & { isExpanded?: boolean }): Color {
        const ret = setInfo<Color>({ type: 'color', defaultValue }, info);
        if (info?.isExpanded) ret.isExpanded = info.isExpanded;
        return ret;
    }

    export interface ColorList extends Base<{ kind: 'interpolate' | 'set', colors: ColorListEntry[] }> {
        type: 'color-list'
        offsets: boolean
        presetKind: 'all' | 'scale' | 'set'
    }
    export function ColorList(defaultValue: { kind: 'interpolate' | 'set', colors: ColorListEntry[] } | ColorListName, info?: Info & { presetKind?: ColorList['presetKind'], offsets?: boolean }): ColorList {
        let def: ColorList['defaultValue'];
        if (typeof defaultValue === 'string') {
            const colors = getColorListFromName(defaultValue);
            def = { kind: colors.type !== 'qualitative' ? 'interpolate' : 'set', colors: colors.list };
        } else {
            def = defaultValue;
        }
        return setInfo<ColorList>({ type: 'color-list', presetKind: info?.presetKind || 'all', defaultValue: def, offsets: !!info?.offsets }, info);
    }

    export interface Vec3 extends Base<Vec3Data>, Range {
        type: 'vec3'
    }
    export function Vec3(defaultValue: Vec3Data, range?: { min?: number, max?: number, step?: number }, info?: Info): Vec3 {
        return setInfo<Vec3>(setRange({ type: 'vec3', defaultValue }, range), info);
    }

    export interface Mat4 extends Base<Mat4Data> {
        type: 'mat4'
    }
    export function Mat4(defaultValue: Mat4Data, info?: Info): Mat4 {
        return setInfo<Mat4>({ type: 'mat4', defaultValue }, info);
    }

    export interface UrlParam extends Base<Asset.Url | string> {
        type: 'url'
    }
    export function Url(url: string | { url: string, body?: string }, info?: Info): UrlParam {
        const defaultValue = typeof url === 'string' ? Asset.Url(url) : Asset.Url(url.url, { body: url.body });
        const ret = setInfo<UrlParam>({ type: 'url', defaultValue }, info);
        return ret;
    }

    export interface FileParam extends Base<Asset.File | null> {
        type: 'file'
        accept?: string
    }
    export function File(info?: Info & { accept?: string, multiple?: boolean }): FileParam {
        const ret = setInfo<FileParam>({ type: 'file', defaultValue: null }, info);
        if (info?.accept) ret.accept = info.accept;
        return ret;
    }

    export interface FileListParam extends Base<Asset.File[] | null> {
        type: 'file-list'
        accept?: string
    }
    export function FileList(info?: Info & { accept?: string, multiple?: boolean }): FileListParam {
        const ret = setInfo<FileListParam>({ type: 'file-list', defaultValue: null }, info);
        if (info?.accept) ret.accept = info.accept;
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
    function setRange<T extends Numeric | Interval | Vec3>(p: T, range?: { min?: number, max?: number, step?: number }) {
        if (!range) return p;
        if (typeof range.min !== 'undefined') p.min = range.min;
        if (typeof range.max !== 'undefined') p.max = range.max;
        if (typeof range.step !== 'undefined') p.step = range.step;
        return p;
    }

    export interface Numeric extends Base<number>, Range {
        type: 'number',
        immediateUpdate?: boolean
    }
    export function Numeric(defaultValue: number, range?: { min?: number, max?: number, step?: number }, info?: Info & { immediateUpdate?: boolean }): Numeric {
        const ret = setInfo<Numeric>(setRange({ type: 'number', defaultValue }, range), info);
        if (info?.immediateUpdate) ret.immediateUpdate = true;
        return ret;
    }

    export interface Interval extends Base<[number, number]>, Range {
        type: 'interval'
    }
    export function Interval(defaultValue: [number, number], range?: { min?: number, max?: number, step?: number }, info?: Info): Interval {
        return setInfo<Interval>(setRange({ type: 'interval', defaultValue }, range), info);
    }

    export interface LineGraph extends Base<Vec2Data[]> {
        type: 'line-graph',
        getVolume?: () => unknown
    }
    export function LineGraph(defaultValue: Vec2Data[], info?: Info & { getVolume?: (binCount?: number) => unknown }): LineGraph {
        const ret = setInfo<LineGraph>({ type: 'line-graph', defaultValue }, info);
        if (info?.getVolume) ret.getVolume = info.getVolume;
        return ret;
    }

    export interface Group<T> extends Base<T> {
        type: 'group',
        params: Params,
        presets?: Select<T>['options'],
        isExpanded?: boolean,
        isFlat?: boolean,
        pivot?: keyof T
    }
    export function Group<T>(params: For<T>, info?: Info & { isExpanded?: boolean, isFlat?: boolean, customDefault?: any, pivot?: keyof T, presets?: Select<T>['options'] }): Group<Normalize<T>> {
        const ret = setInfo<Group<Normalize<T>>>({ type: 'group', defaultValue: info?.customDefault || getDefaultValues(params as any as Params) as any, params: params as any as Params }, info);
        if (info?.presets) ret.presets = info.presets;
        if (info?.isExpanded) ret.isExpanded = info.isExpanded;
        if (info?.isFlat) ret.isFlat = info.isFlat;
        if (info?.pivot) ret.pivot = info.pivot as any;
        return ret;
    }
    export function EmptyGroup(info?: Info) {
        return Group({}, info);
    }

    export interface NamedParams<T = any, K = string> { name: K, params: T }
    export type NamedParamUnion<P extends Params, K extends keyof P = keyof P> = K extends any ? NamedParams<P[K]['defaultValue'], K> : never
    export interface Mapped<T extends NamedParams<any, any>> extends Base<T> {
        type: 'mapped',
        select: Select<string>,
        map(name: string): Any
    }
    export function Mapped<T>(defaultKey: string, names: ([string, string] | [string, string, string])[], map: (name: string) => Any, info?: Info & { cycle?: boolean }): Mapped<NamedParams<T>> {
        const name = checkDefaultKey(defaultKey, names);
        return setInfo<Mapped<NamedParams<T>>>({
            type: 'mapped',
            defaultValue: { name, params: map(name).defaultValue as any },
            select: Select<string>(name, names, info),
            map
        }, info);
    }
    export function MappedStatic<C extends Params>(defaultKey: keyof C, map: C, info?: Info & { options?: [keyof C, string][], cycle?: boolean }): Mapped<NamedParamUnion<C>> {
        const options: [string, string][] = info?.options
            ? info.options as [string, string][]
            : Object.keys(map).map(k => [k, map[k].label || stringToWords(k)]) as [string, string][];
        const name = checkDefaultKey(defaultKey, options);
        return setInfo<Mapped<NamedParamUnion<C>>>({
            type: 'mapped',
            defaultValue: { name, params: map[name].defaultValue } as any,
            select: Select<string>(name as string, options, info),
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
        return setInfo<ObjectList<Normalize<T>>>({ type: 'object-list', element: element as any as Params, getLabel, ctor: _defaultObjectListCtor, defaultValue: (info?.defaultValue) || [] }, info);
    }
    function _defaultObjectListCtor(this: ObjectList) { return getDefaultValues(this.element) as any; }


    function unsetGetValue() {
        throw new Error('getValue not set. Fix runtime.');
    }

    // getValue needs to be assigned by a runtime because it might not be serializable
    export interface ValueRef<T = any> extends Base<{ ref: string, getValue: () => T }> {
        type: 'value-ref',
        resolveRef: (ref: string, getData: (ref: string) => any) => T,
        // a provider because the list changes over time
        getOptions: (ctx: any) => Select<string>['options'],
    }
    export function ValueRef<T>(getOptions: ValueRef['getOptions'], resolveRef: ValueRef<T>['resolveRef'], info?: Info & { defaultRef?: string }) {
        return setInfo<ValueRef<T>>({ type: 'value-ref', defaultValue: { ref: info?.defaultRef ?? '', getValue: unsetGetValue as any }, getOptions, resolveRef }, info);
    }

    export interface DataRef<T = any> extends Base<{ ref: string, getValue: () => T }> {
        type: 'data-ref'
    }
    export function DataRef<T>(info?: Info & { defaultRef?: string }) {
        return setInfo<DataRef<T>>({ type: 'data-ref', defaultValue: { ref: info?.defaultRef ?? '', getValue: unsetGetValue as any } }, info);
    }

    export interface Converted<T, C> extends Base<T> {
        type: 'converted',
        converted: Any,
        /** converts from prop value to display value */
        fromValue(v: T): C,
        /** converts from display value to prop value */
        toValue(v: C): T
    }
    export function Converted<T, C extends Any>(fromValue: (v: T) => C['defaultValue'], toValue: (v: C['defaultValue']) => T, converted: C): Converted<T, C['defaultValue']> {
        return setInfo({ type: 'converted', defaultValue: toValue(converted.defaultValue), converted, fromValue, toValue }, converted);
    }

    export interface Conditioned<T, P extends Base<T>, C = { [k: string]: P }> extends Base<T> {
        type: 'conditioned',
        select: Select<string>,
        conditionParams: C
        conditionForValue(v: T): keyof C
        conditionedValue(v: T, condition: keyof C): T,
    }
    export function Conditioned<T, P extends Base<T>, C extends {} = { [k: string]: P }>(defaultValue: T, conditionParams: C, conditionForValue: (v: T) => keyof C, conditionedValue: (v: T, condition: keyof C) => T, info?: Info): Conditioned<T, P, C> {
        const options = Object.keys(conditionParams).map(k => [k, k]) as [string, string][];
        return setInfo({ type: 'conditioned', select: Select<string>(conditionForValue(defaultValue) as string, options, info), defaultValue, conditionParams, conditionForValue, conditionedValue }, info);
    }

    export interface Script extends Base<ScriptData> {
        type: 'script'
    }
    export function Script(defaultValue: Script['defaultValue'], info?: Info): Script {
        return setInfo<Script>({ type: 'script', defaultValue }, info);
    }

    export type Any =
        | Value<any> | Select<any> | MultiSelect<any> | BooleanParam | Text | Color | Vec3 | Mat4 | Numeric | FileParam | UrlParam | FileListParam | Interval | LineGraph
        | ColorList | Group<any> | Mapped<any> | Converted<any, any> | Conditioned<any, any, any> | Script | ObjectList | ValueRef | DataRef

    export type Params = { [k: string]: Any }
    export type Values<T extends Params = Params> = { [k in keyof T]: T[k]['defaultValue'] }
    /** This is required for params with optional values */
    export type ValuesFor<T extends For<any>> = Normalize<{ [k in keyof T]: T[k]['defaultValue'] }>

    type Optionals<P> = { [K in keyof P]-?: undefined extends P[K] ? K : never }[keyof P]
    type NonOptionals<P> = { [K in keyof P]-?: undefined extends P[K] ? never : K }[keyof P]
    export type Normalize<P> = Pick<P, NonOptionals<P>> & Partial<Pick<P, Optionals<P>>>
    export type For<P> = { [K in keyof P]-?: Base<P[K]> }
    export type Def<P> = { [K in keyof P]: Any }


    export function For<P>(params: For<P>): For<P> {
        return 0 as any;
    }

    export function getDefaultValues<T extends Params>(params: T) {
        const d: { [k: string]: any } = {};
        for (const k of Object.keys(params)) {
            if (params[k].isOptional) continue;
            d[k] = params[k].defaultValue;
        }
        return d as Values<T>;
    }

    function _resolveRef(resolve: (ref: string, getData: (ref: string) => any) => any, ref: string, getData: (ref: string) => any) {
        return () => resolve(ref, getData);
    }

    function resolveRefValue(p: Any, value: any, getData: (ref: string) => any) {
        if (!value) return;

        if (p.type === 'value-ref') {
            const v = value as ValueRef['defaultValue'];
            if (!v.ref) v.getValue = () => { throw new Error('Unset ref in ValueRef value.'); };
            else v.getValue = _resolveRef(p.resolveRef, v.ref, getData);
        } else if (p.type === 'data-ref') {
            const v = value as ValueRef['defaultValue'];
            if (!v.ref) v.getValue = () => { throw new Error('Unset ref in ValueRef value.'); };
            else v.getValue = _resolveRef(getData, v.ref, getData);
        } else if (p.type === 'group') {
            resolveRefs(p.params, value, getData);
        } else if (p.type === 'mapped') {
            const v = value as NamedParams;
            const param = p.map(v.name);
            resolveRefValue(param, v.params, getData);
        } else if (p.type === 'object-list') {
            if (!hasValueRef(p.element)) return;
            for (const e of value) {
                resolveRefs(p.element, e, getData);
            }
        }
    }

    function hasParamValueRef(p: Any) {
        if (p.type === 'value-ref' || p.type === 'data-ref') {
            return true;
        } else if (p.type === 'group') {
            if (hasValueRef(p.params)) return true;
        } else if (p.type === 'mapped') {
            for (const [o] of p.select.options) {
                if (hasParamValueRef(p.map(o))) return true;
            }
        } else if (p.type === 'object-list') {
            return hasValueRef(p.element);
        }
        return false;
    }

    function hasValueRef(params: Params) {
        for (const n of Object.keys(params)) {
            if (hasParamValueRef(params[n])) return true;
        }
        return false;
    }

    export function resolveRefs(params: Params, values: any, getData: (ref: string) => any) {
        for (const n of Object.keys(params)) {
            resolveRefValue(params[n], values?.[n], getData);
        }
    }

    export function setDefaultValues<T extends Params>(params: T, defaultValues: Values<T>) {
        for (const k of Object.keys(params)) {
            if (params[k].isOptional) continue;
            params[k].defaultValue = defaultValues[k];
        }
    }

    export function clone<P extends Params>(params: P): P {
        return deepClone(params);
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

        if (typeof a !== 'object' || typeof b !== 'object') return false;
        for (const k of Object.keys(params)) {
            if (!isParamEqual(params[k], a[k], b[k])) return false;
        }
        return true;
    }

    export function isParamEqual(p: Any, a: any, b: any): boolean {
        if (a === b) return true;

        if (p.type === 'group') {
            return areEqual(p.params, a, b);
        } else if (p.type === 'mapped') {
            const u = a as NamedParams, v = b as NamedParams;
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
        } else if (p.type === 'mat4') {
            return Mat4Data.areEqual(a, b, EPSILON);
        } else if (p.type === 'script') {
            const u = a as Script['defaultValue'], v = b as Script['defaultValue'];
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
            return shallowEqualObjects(a, b);
        }

        // a === b was checked at the top.
        return false;
    }

    export function merge<P extends Params>(params: P, a: any, b: any): Values<P> {
        if (a === undefined) return { ...b };
        if (b === undefined) return { ...a };

        const o = Object.create(null);
        for (const k of Object.keys(params)) {
            o[k] = mergeParam(params[k], a[k], b[k]);
        }
        return o;
    }

    export function mergeParam(p: Any, a: any, b: any): any {
        if (a === undefined) return typeof b === 'object' && !Array.isArray(b) ? { ...b } : b;
        if (b === undefined) return typeof a === 'object' && !Array.isArray(a) ? { ...a } : a;

        if (p.type === 'group') {
            return merge(p.params, a, b);
        } else if (p.type === 'mapped') {
            const u = a as NamedParams, v = b as NamedParams;
            if (u.name !== v.name) return { ...v };
            const map = p.map(v.name);
            return {
                name: v.name,
                params: mergeParam(map, u.params, v.params)
            };
        } else if (p.type === 'value') {
            return b;
        } else if (typeof a === 'object' && typeof b === 'object') {
            if (Array.isArray(b)) {
                return b;
            }
            return { ...a, ...b };
        } else {
            return b;
        }
    }

    function selectHasOption(p: Select<any> | MultiSelect<any>, v: any) {
        for (const o of p.options) {
            if (o[0] === v) return true;
        }
        return false;
    }

    function normalizeParam(p: Any, value: any, defaultIfUndefined: boolean): any {
        if (value === void 0 || value === null) {
            return defaultIfUndefined ? p.defaultValue : void 0;
        }

        // TODO: is this a good idea and will work well?
        // if (typeof p.defaultValue !== typeof value) {
        //     return p.defaultValue;
        // }

        if (p.type === 'value') {
            return value;
        } else if (p.type === 'group') {
            const ret = Object.create(null);
            for (const key of Object.keys(p.params)) {
                const param = p.params[key];
                if (value[key] === void 0) {
                    if (defaultIfUndefined) ret[key] = param.defaultValue;
                } else {
                    ret[key] = normalizeParam(param, value[key], defaultIfUndefined);
                }
            }
            return ret;
        } else if (p.type === 'mapped') {
            const v = value as NamedParams;
            if (typeof v.name !== 'string') {
                return p.defaultValue;
            }
            if (typeof v.params === 'undefined') {
                return defaultIfUndefined ? p.defaultValue : void 0;
            }

            if (!selectHasOption(p.select, v.name)) {
                return p.defaultValue;
            }

            const param = p.map(v.name);
            return {
                name: v.name,
                params: normalizeParam(param, v.params, defaultIfUndefined)
            };
        } else if (p.type === 'select') {
            if (!selectHasOption(p, value)) return p.defaultValue;
            return value;
        } else if (p.type === 'multi-select') {
            if (!Array.isArray(value)) return p.defaultValue;
            const ret = value.filter(function (this: MultiSelect<any>, v: any) { return selectHasOption(this, v); }, p);
            if (value.length > 0 && ret.length === 0) return p.defaultValue;
            return ret;
        } else if (p.type === 'object-list') {
            if (!Array.isArray(value)) return p.defaultValue;
            return value.map(v => normalizeParams(p.element, v, defaultIfUndefined ? 'all' : 'skip'));
        }

        // TODO: validate/normalize all param types "properly"??

        return value;
    }

    export function normalizeParams(p: Params, value: any, defaultIfUndefined: 'all' | 'children' | 'skip') {
        if (typeof value !== 'object' || value === null) {
            return defaultIfUndefined ? getDefaultValues(p) : value;
        }

        const ret = Object.create(null);
        for (const key of Object.keys(p)) {
            const param = p[key];
            if (value[key] === void 0) {
                if (defaultIfUndefined === 'all') ret[key] = param.defaultValue;
            } else {
                ret[key] = normalizeParam(param, value[key], defaultIfUndefined !== 'skip');
            }
        }
        return ret;
    }

    /**
     * Map an object to a list of [K, string][] to be used as options, stringToWords for key used by default (or identity of null).
     *
     * if options is { [string]: string } and mapping is not provided, use the Value.
     */
    export function objectToOptions<K extends string, V>(options: { [k in K]: V }, f?: null | ((k: K, v: V) => string | [string, string])): [K, string][] {
        const ret: ([K, string] | [K, string, string])[] = [];
        for (const k of Object.keys(options) as K[]) {
            if (!f) {
                if (typeof options[k as K] === 'string') ret.push([k as K, options[k as K] as any]);
                else ret.push([k as K, f === null ? k : stringToWords(k)]);
            } else {
                const o = f(k as K, options[k as K]);
                ret.push(typeof o === 'string' ? [k, o] : [k, o[0], o[1]]);
            }
        }
        return ret as [K, string][];
    }

    /**
     * Map array of options using stringToWords by default (or identity of null).
     */
    export function arrayToOptions<V extends string>(xs: readonly V[], f?: null | ((v: V) => string)): [V, string][] {
        const ret: [V, string][] = [];
        for (const x of xs) {
            if (!f) {
                ret.push([x, f === null ? x : stringToWords(x)]);
            } else {
                ret.push([x, f(x)]);
            }
        }
        return ret;
    }

    export function optionLabel<T>(param: Select<T>, value: T) {
        for (const o of param.options) {
            if (o[0] === value) return o[1];
        }
        return '';
    }

    function checkDefaultKey<T>(k: T, options: readonly (readonly [T, string] | readonly [T, string, string | undefined])[]) {
        for (const o of options) {
            if (o[0] === k) return k;
        }
        return options.length > 0 ? options[0][0] : void 0 as any as T;
    }

    export function withDefaults<T extends Params>(schema: T, updates: Partial<ValuesFor<T>>): T {
        const next: any = {};
        for (const k of Object.keys(updates)) {
            const v = (updates as any)[k];
            if (!schema[k] || v === null || v === undefined || schema[k].defaultValue === v) continue;
            next[k] = { ...schema[k], defaultValue: v };
        }
        return { ...schema, ...next };
    }
}