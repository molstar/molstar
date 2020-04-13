/*
 * Copyright (c) 2018 Mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

type Type<T = any> =
    | Type.Any | Type.AnyValue | Type.Variable<T> | Type.Value<T>
    | Type.Container<T> | Type.Union<T> | Type.OneOf<T>

namespace Type {
    export interface Any { kind: 'any', '@type': any }
    export interface Variable<T> { kind: 'variable', name: string, type: Type, isConstraint: boolean, '@type': any }
    export interface AnyValue { kind: 'any-value', '@type': any }
    export interface Value<T> { kind: 'value', namespace: string, name: string, parent?: Value<any>, '@type': T }
    export interface Container<T> { kind: 'container', namespace: string, name: string, alias?: string, child: Type, '@type': T }
    export interface Union<T> { kind: 'union', types: Type[], '@type': T }
    export interface OneOf<T> { kind: 'oneof', type: Value<T>, namespace: string, name: string, values: { [v: string]: boolean | undefined }, '@type': T }

    export function Variable<T = any>(name: string, type: Type, isConstraint?: boolean): Variable<T> { return { kind: 'variable', name, type: type, isConstraint } as any; }
    export function Value<T>(namespace: string, name: string, parent?: Value<any>): Value<T> { return { kind: 'value', namespace, name, parent } as any; }
    export function Container<T = any>(namespace: string, name: string, child: Type, alias?: string): Container<T> { return { kind: 'container', namespace, name, child, alias } as any; }
    export function Union<T = any>(types: Type[]): Union<T> { return { kind: 'union', types } as any; }
    export function OneOf<T = any>(namespace: string, name: string, type: Value<T>, values: any[]): OneOf<T> {
        const vals = Object.create(null);
        for (const v of values) vals[v] = true;
        return { kind: 'oneof', namespace, name, type, values: vals } as any;
    }

    export const Any: Any = { kind: 'any' } as any;
    export const AnyValue: AnyValue = { kind: 'any-value' } as any;

    export const Num = Value<number>('', 'Number');
    export const Str = Value<string>('', 'String');
    export const Bool = OneOf<boolean>('', 'Bool', Str as any, ['true', 'false']);

    export function oneOfValues({ values }: OneOf<any>) { return Object.keys(values).sort(); }
}

export default Type;