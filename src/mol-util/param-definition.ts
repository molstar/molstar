/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Color as ColorData } from './color';
import { shallowClone } from 'mol-util';

export namespace ParamDefinition {
    export interface Base<T> {
        label: string
        description: string
        defaultValue: T
    }

    export interface GenericValue<T> extends Base<T> {
        type: 'generic-value'
    }
    export function GenericValue<T>(label: string, description: string, defaultValue: T): GenericValue<T> {
        return { type: 'generic-value', label, description, defaultValue }
    }

    export interface Select<T extends string> extends Base<T> {
        type: 'select'
        /** array of (value, label) tuples */
        options: [T, string][]
    }
    export function Select<T extends string>(label: string, description: string, defaultValue: T, options: [T, string][]): Select<T> {
        return { type: 'select', label, description, defaultValue, options }
    }

    export interface MultiSelect<E extends string, T = E[]> extends Base<T> {
        type: 'multi-select'
        /** array of (value, label) tuples */
        options: [E, string][]
    }
    export function MultiSelect<E extends string, T = E[]>(label: string, description: string, defaultValue: T, options: [E, string][]): MultiSelect<E, T> {
        return { type: 'multi-select', label, description, defaultValue, options }
    }

    export interface Boolean extends Base<boolean> {
        type: 'boolean'
    }
    export function Boolean(label: string, description: string, defaultValue: boolean): Boolean {
        return { type: 'boolean', label, description, defaultValue }
    }

    export interface Range extends Base<number> {
        type: 'range'
        min: number
        max: number
        /** if an `integer` parse value with parseInt, otherwise use parseFloat */
        step: number
    }
    export function Range(label: string, description: string, defaultValue: number, min: number, max: number, step: number): Range {
        return { type: 'range', label, description, defaultValue, min, max, step }
    }

    export interface Text extends Base<string> {
        type: 'text'
    }
    export function Text(label: string, description: string, defaultValue: string = ''): Text {
        return { type: 'text', label, description, defaultValue }
    }

    export interface Color extends Base<ColorData> {
        type: 'color'
    }
    export function Color(label: string, description: string, defaultValue: ColorData): Color {
        return { type: 'color', label, description, defaultValue }
    }

    export interface Numeric extends Base<number> {
        type: 'number'
        min: number
        max: number
        /** if an `integer` parse value with parseInt, otherwise use parseFloat */
        step: number
    }
    export function Numeric(label: string, description: string, defaultValue: number, min: number, max: number, step: number): Numeric {
        return { type: 'number', label, description, defaultValue, min, max, step }
    }

    export interface Interval extends Base<[number, number]> {
        type: 'interval'
    }
    export function Interval(label: string, description: string, defaultValue: [number, number]): Interval {
        return { type: 'interval', label, description, defaultValue }
    }

    export interface Group<T extends { [key: string]: any }> extends Base<T> {
        type: 'group',
        params: T
    }
    export function Group<T extends { [key: string]: any }>(label: string, description: string, params: T): Group<T> {
        return { type: 'group', label, description, defaultValue: getDefaultValues(params), params };
    }

    export interface Mapped<T> extends Base<{ name: string, params: T }> {
        type: 'mapped',
        select: Select<string>,
        map(name: string): Params
    }
    export function Mapped<T>(label: string, description: string, defaultKey: string, names: [string, string][], map: Mapped<T>['map']): Mapped<T> {
        return {
            type: 'mapped',
            label,
            description,
            defaultValue: { name: defaultKey, params: getDefaultValues(map(defaultKey)) as any },
            select: Select<string>(label, description, defaultKey, names),
            map
        };
    }

    export type Any = Select<any> | MultiSelect<any> | Boolean | Range | Text | Color | Numeric | Interval | Group<any> | Mapped<any>

    export type Params = { [k: string]: Any }
    export type Values<T extends Params> = { [k in keyof T]: T[k]['defaultValue'] }

    export function getDefaultValues<T extends Params>(params: T) {
        const d: { [k: string]: any } = {}
        Object.keys(params).forEach(k => d[k] = params[k].defaultValue)
        return d as Values<T>
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