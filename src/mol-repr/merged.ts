/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { GraphicsRenderObject, MergedGraphicsRenderObject, canMergeRenderObjects, createMergedRenderObject } from '../mol-gl/render-object';
import { RenderMergedValue, resolveRenderMerged } from '../mol-geo/geometry/base';
import { PickingId } from '../mol-geo/geometry/picking';
import { Loci, EmptyLoci } from '../mol-model/loci';

/** One visual and its renderObject, ready to be grouped/merged by `mergeRenderObjectsByType`. */
export type MergeableEntry<V> = { visual: V, renderObject: GraphicsRenderObject }

/** A currently-merged group of render objects, keyed by their shared geometry type. */
export type MergedEntry<V> = { object: MergedGraphicsRenderObject, visuals: V[], key: string }

/**
 * Groups `entries` by `renderObject.type` and merges each group of >= 2 mergeable
 * render objects into one shared render item, when `renderMerged` resolves to true
 * for the group's size. Reuses `prevMergedByType`'s entry for a group whose member
 * ids are unchanged, instead of rebuilding a merged render object every update.
 *
 * Clears and refills `renderObjects` with the merged and any left-over unmerged
 * render objects. Returns the new by-type map, to pass as `prevMergedByType` next time.
 */
export function mergeRenderObjectsByType<V>(
    entries: readonly MergeableEntry<V>[],
    renderMerged: RenderMergedValue,
    prevMergedByType: ReadonlyMap<string, MergedEntry<V>>,
    renderObjects: GraphicsRenderObject[]
): Map<string, MergedEntry<V>> {
    renderObjects.length = 0;
    const nextMerged = new Map<string, MergedEntry<V>>();
    const consumed = new Set<GraphicsRenderObject>();

    if (renderMerged !== false) {
        const byType = new Map<string, MergeableEntry<V>[]>();
        for (const e of entries) {
            const group = byType.get(e.renderObject.type);
            if (group) group.push(e);
            else byType.set(e.renderObject.type, [e]);
        }
        byType.forEach((group, type) => {
            if (group.length < 2 || !resolveRenderMerged(renderMerged, group.length)) return;
            if (!canMergeRenderObjects(group.map(e => e.renderObject))) return;

            const key = group.map(e => e.renderObject.id).join(',');
            const prev = prevMergedByType.get(type);
            const entry = prev && prev.key === key ? prev : {
                object: createMergedRenderObject(group.map(e => e.renderObject)),
                visuals: group.map(e => e.visual),
                key,
            };
            nextMerged.set(type, entry);
            renderObjects.push(entry.object);
            for (const e of group) consumed.add(e.renderObject);
        });
    }

    for (const e of entries) {
        if (!consumed.has(e.renderObject)) renderObjects.push(e.renderObject);
    }

    return nextMerged;
}

/**
 * Resolves a `pickingId` that may target a merged render object: finds the
 * `MergedEntry` whose object id matches, translates the global instance id to
 * the owning member's index and local instance id, and delegates to `resolve`.
 * Returns undefined when `pickingId` does not target any currently merged
 * object, letting the caller fall back to per-visual picking.
 */
export function pickMergedLoci<V>(
    mergedByType: ReadonlyMap<string, MergedEntry<V>>,
    pickingId: PickingId,
    resolve: (visual: V, index: number, localInstanceId: number, object: MergedGraphicsRenderObject) => Loci
): Loci | undefined {
    let loci: Loci | undefined;
    mergedByType.forEach(({ object, visuals }) => {
        if (loci !== undefined || pickingId.objectId !== object.id) return;
        const { segments } = object.merged;
        const { instanceId } = pickingId;
        for (let i = 0, il = segments.count; i < il; ++i) {
            const base = segments.instanceBases[i];
            if (instanceId >= base && instanceId < base + segments.instanceCounts[i]) {
                loci = resolve(visuals[i], i, instanceId - base, object);
                return;
            }
        }
        loci = EmptyLoci;
    });
    return loci;
}
