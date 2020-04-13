/*
 * Copyright (c) 2018 Mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import Type from './type';
import { MSymbol, Arguments, isSymbol } from './symbol';

export function symbol<A extends Arguments, T extends Type<S>, S>(args: A, type: T, description?: string) {
    return MSymbol('', args, type, description);
}

export function normalizeTable(table: any) {
    _normalizeTable('', '', table);
}

export function symbolList(table: any): MSymbol[] {
    const list: MSymbol[] = [];
    _symbolList(table, list);
    return list;
}

function formatKey(key: string) {
    const regex = /([a-z])([A-Z])([a-z]|$)/g;
    // do this twice because 'xXxX'
    return key.replace(regex, (s, a, b, c) => `${a}-${b.toLocaleLowerCase()}${c}`).replace(regex, (s, a, b, c) => `${a}-${b.toLocaleLowerCase()}${c}`);
}

function _normalizeTable(namespace: string, key: string, obj: any) {
    if (isSymbol(obj)) {
        obj.info.namespace = namespace;
        obj.info.name = obj.info.name || formatKey(key);
        obj.id = `${obj.info.namespace}.${obj.info.name}`;
        return;
    }
    const currentNs = `${obj['@namespace'] || formatKey(key)}`;
    const newNs = namespace ? `${namespace}.${currentNs}` : currentNs;
    for (const childKey of Object.keys(obj)) {
        if (typeof obj[childKey] !== 'object' && !isSymbol(obj[childKey])) continue;
        _normalizeTable(newNs, childKey, obj[childKey]);
    }
}

function _symbolList(obj: any, list: MSymbol[]) {
    if (isSymbol(obj)) {
        list.push(obj);
        return;
    }
    for (const childKey of Object.keys(obj)) {
        if (typeof obj[childKey] !== 'object' && !isSymbol(obj[childKey])) continue;
        _symbolList(obj[childKey], list);
    }
}