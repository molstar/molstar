/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 * 
 * @author Taylor Hoffmann <taylor@hoffmann.io>
 */

import { Structure } from '../../mol-model/structure';
import { Vec3 } from '../../mol-math/linear-algebra';

const serialMapCache = new WeakMap<Structure, Map<string, Map<string, number>>>();
const operatorHklCache = new WeakMap<Structure, { min: Vec3, max: Vec3, map: Map<string, number> }>();

export function getCachedSerialMap(structure: Structure, key: string, factory: () => Map<string, number>): Map<string, number> {
    let byKey = serialMapCache.get(structure);
    if (!byKey) {
        byKey = new Map();
        serialMapCache.set(structure, byKey);
    }
    let map = byKey.get(key);
    if (!map) {
        map = factory();
        byKey.set(key, map);
    }
    return map;
}

export function getCachedOperatorHklMap(structure: Structure, factory: () => { min: Vec3, max: Vec3, map: Map<string, number> }) {
    let data = operatorHklCache.get(structure);
    if (!data) {
        data = factory();
        operatorHklCache.set(structure, data);
    }
    return data;
}
