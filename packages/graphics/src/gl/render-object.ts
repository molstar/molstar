/**
 * Copyright (c) 2018-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import type { RenderableState, Renderable } from './renderable.js';
import { idFactory } from '@molstar/core/util/id-factory';
import type { WebGLContext } from './webgl/context.js';
import { type DirectVolumeValues, DirectVolumeRenderable } from './renderable/direct-volume.js';
import { type MeshValues, MeshRenderable } from './renderable/mesh.js';
import { type PointsValues, PointsRenderable } from './renderable/points.js';
import { type LinesValues, LinesRenderable } from './renderable/lines.js';
import { type SpheresValues, SpheresRenderable } from './renderable/spheres.js';
import { type TextValues, TextRenderable } from './renderable/text.js';
import { type TextureMeshValues, TextureMeshRenderable } from './renderable/texture-mesh.js';
import { type ImageValues, ImageRenderable } from './renderable/image.js';
import { CylindersRenderable, type CylindersValues } from './renderable/cylinders.js';
import type { Transparency } from './webgl/render-item.js';
import type { GlobalDefines } from './renderable/schema.js';
import { assertUnreachable } from '@molstar/core/util/type-helpers';
import { canMergeValues, createMergedValues, isMergeableType, type Merged, type MergeableValues, MergedRenderable } from './renderable/merged.js';

const getNextId = idFactory(0, 0x7FFFFFFF);

export const getNextMaterialId = idFactory(0, 0x7FFFFFFF);

export interface GraphicsRenderObject<T extends RenderObjectType = RenderObjectType> {
    readonly id: number,
    readonly type: T,
    readonly values: RenderObjectValues<T>,
    readonly state: RenderableState,
    readonly materialId: number
}

export type RenderObjectType = 'mesh' | 'points' | 'spheres' | 'cylinders' | 'text' | 'lines' | 'direct-volume' | 'image' | 'texture-mesh'

export type RenderObjectValues<T extends RenderObjectType> =
    T extends 'mesh' ? MeshValues :
        T extends 'points' ? PointsValues :
            T extends 'spheres' ? SpheresValues :
                T extends 'cylinders' ? CylindersValues :
                    T extends 'text' ? TextValues :
                        T extends 'lines' ? LinesValues :
                            T extends 'direct-volume' ? DirectVolumeValues :
                                T extends 'image' ? ImageValues :
                                    T extends 'texture-mesh' ? TextureMeshValues : never

//

export function createRenderObject<T extends RenderObjectType>(type: T, values: RenderObjectValues<T>, state: RenderableState, materialId: number): GraphicsRenderObject<T> {
    return { id: getNextId(), type, values, state, materialId } as GraphicsRenderObject<T>;
}

//

/**
 * A render object that shares a single render item between the values of
 * multiple member render objects to reduce draw calls. The members are
 * not meant to be added to a scene themselves.
 */
export interface MergedGraphicsRenderObject<T extends RenderObjectType = RenderObjectType> extends GraphicsRenderObject<T> {
    readonly members: ReadonlyArray<GraphicsRenderObject<T>>
    readonly merged: Merged
}

export function isMergedRenderObject<T extends RenderObjectType>(o: GraphicsRenderObject<T>): o is MergedGraphicsRenderObject<T> {
    return 'members' in o;
}

export function canMergeRenderObjects(members: readonly GraphicsRenderObject[]): boolean {
    if (members.length < 2) return false;
    const { type, materialId } = members[0];
    if (!isMergeableType(type)) return false;
    for (const m of members) {
        if (m.type !== type) return false;
        if (m.materialId !== materialId) return false;
    }
    if (!canMergeValues(type, members.map(m => m.values as MergeableValues))) return false;
    return true;
}

/**
 * Creates a merged render object from member render objects, assumes
 * `canMergeRenderObjects(members)`. Shares the state of the first member.
 */
export function createMergedRenderObject<T extends RenderObjectType>(members: readonly GraphicsRenderObject<T>[]): MergedGraphicsRenderObject<T> {
    if (!isMergeableType(members[0].type)) throw new Error(`unsupported merged render object type '${members[0].type}'`);
    const merged = createMergedValues(members[0].type, members.map(m => m.values as MergeableValues));
    return { id: getNextId(), type: members[0].type, values: merged.values as unknown as RenderObjectValues<T>, state: members[0].state, materialId: members[0].materialId, members, merged };
}

export function createRenderable<T extends RenderObjectType>(ctx: WebGLContext, o: GraphicsRenderObject<T>, transparency: Transparency, globals: GlobalDefines): Renderable<any> {
    if (isMergedRenderObject(o)) {
        return MergedRenderable(ctx, o.id, o.merged, o.members.map(m => m.values as MergeableValues), o.state, o.materialId, transparency, globals);
    }
    switch (o.type) {
        case 'mesh': return MeshRenderable(ctx, o.id, o.values as MeshValues, o.state, o.materialId, transparency, globals);
        case 'points': return PointsRenderable(ctx, o.id, o.values as PointsValues, o.state, o.materialId, transparency, globals);
        case 'spheres': return SpheresRenderable(ctx, o.id, o.values as SpheresValues, o.state, o.materialId, transparency, globals);
        case 'cylinders': return CylindersRenderable(ctx, o.id, o.values as CylindersValues, o.state, o.materialId, transparency, globals);
        case 'text': return TextRenderable(ctx, o.id, o.values as TextValues, o.state, o.materialId, transparency, globals);
        case 'lines': return LinesRenderable(ctx, o.id, o.values as LinesValues, o.state, o.materialId, transparency, globals);
        case 'direct-volume': return DirectVolumeRenderable(ctx, o.id, o.values as DirectVolumeValues, o.state, o.materialId, transparency, globals);
        case 'image': return ImageRenderable(ctx, o.id, o.values as ImageValues, o.state, o.materialId, transparency, globals);
        case 'texture-mesh': return TextureMeshRenderable(ctx, o.id, o.values as TextureMeshValues, o.state, o.materialId, transparency, globals);
    }
    assertUnreachable(o.type);
}
