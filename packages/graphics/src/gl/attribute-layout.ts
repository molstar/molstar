/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import type { AttributeKind, AttributeItemSize } from './buffer.js';
import type { AttributeSpec, RenderableSchema } from './renderable/schema.js';

/**
 * Backend-neutral attribute locations of a schema. Locations are fixed per schema, so they
 * are the same for all programs (variants) of a render item. They become explicit
 * `@location`s / vertex buffer layouts on WebGPU.
 */

export type AttributeLayoutEntry = {
    readonly name: string;
    readonly kind: AttributeKind;
    readonly itemSize: AttributeItemSize;
    readonly divisor: number;
    readonly location: number;
    /** 4 for `mat4` (one per column), 1 otherwise */
    readonly locationCount: number;
};

export interface AttributeLayout {
    readonly entries: readonly AttributeLayoutEntry[];
    readonly byName: { readonly [k: string]: AttributeLayoutEntry };
    /** number of locations used */
    readonly size: number;
    readonly hash: string;
}

/**
 * Locations in schema order, except that location 0 goes to the first per-vertex attribute
 * without `usedIf` (or else to the first attribute without `usedIf`): some drivers emulate
 * a disabled or instanced attribute 0.
 */
export function getAttributeLayout(schema: RenderableSchema): AttributeLayout {
    // assumes attribute specs are not added to a schema after first use
    let layout = layoutCache.get(schema);
    if (layout === undefined) {
        layout = computeAttributeLayout(schema);
        layoutCache.set(schema, layout);
    }
    return layout;
}

const layoutCache = new WeakMap<RenderableSchema, AttributeLayout>();

function computeAttributeLayout(schema: RenderableSchema): AttributeLayout {
    const specs: [string, AttributeSpec<AttributeKind>][] = [];
    for (const k of Object.keys(schema)) {
        const spec = schema[k];
        if (spec.type === 'attribute') specs.push([k, spec]);
    }
    let index = specs.findIndex(([, s]) => !s.usedIf && s.divisor === 0);
    if (index === -1) index = specs.findIndex(([, s]) => !s.usedIf);
    if (index > 0) specs.unshift(...specs.splice(index, 1));

    const entries: AttributeLayoutEntry[] = [];
    const byName: { [k: string]: AttributeLayoutEntry } = {};
    let location = 0;
    for (const [name, spec] of specs) {
        const locationCount = spec.itemSize === 16 ? 4 : 1;
        const entry = { name, kind: spec.kind, itemSize: spec.itemSize, divisor: spec.divisor, location, locationCount };
        entries.push(entry);
        byName[name] = entry;
        location += locationCount;
    }
    return {
        entries,
        byName,
        size: location,
        hash: entries.map((e) => `${e.name}:${e.location}`).join(','),
    };
}
