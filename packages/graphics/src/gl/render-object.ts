/**
 * Copyright (c) 2018-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import type { RenderableState, Renderable } from './renderable.js';
import { idFactory } from '@molstar/core/util/id-factory';
import type { GraphicsContext } from './context.js';
import {
  type DirectVolumeValues,
  DirectVolumeRenderable,
  DirectVolumeRenderableSchema,
} from './renderable/direct-volume.js';
import { type MeshValues, MeshRenderable, MeshRenderableSchema } from './renderable/mesh.js';
import { type PointsValues, PointsRenderable, PointsRenderableSchema } from './renderable/points.js';
import { type LinesValues, LinesRenderable, LinesRenderableSchema } from './renderable/lines.js';
import { type SpheresValues, SpheresRenderable, SpheresRenderableSchema } from './renderable/spheres.js';
import { type TextValues, TextRenderable, TextRenderableSchema } from './renderable/text.js';
import {
  type TextureMeshValues,
  TextureMeshRenderable,
  TextureMeshRenderableSchema,
} from './renderable/texture-mesh.js';
import { type ImageValues, ImageRenderable, ImageRenderableSchema } from './renderable/image.js';
import { CylindersRenderable, CylindersRenderableSchema, type CylindersValues } from './renderable/cylinders.js';
import type { Transparency } from './render-item.js';
import type { GlobalDefines, RenderableSchema } from './renderable/schema.js';
import {
  CylindersShaderCode,
  DirectVolumeShaderCode,
  ImageShaderCode,
  LinesShaderCode,
  MeshShaderCode,
  PointsShaderCode,
  type ShaderCode,
  SpheresShaderCode,
  TextShaderCode,
} from './shader-code.js';
import { assertUnreachable } from '@molstar/core/util/type-helpers';
import { isNullTexture, type Texture, transferTexture } from './texture.js';
import { ValueCell } from '@molstar/core/util/value-cell';
import {
  canMergeValues,
  createMergedValues,
  isMergeableType,
  type Merged,
  type MergeableValues,
  MergedRenderable,
  type MergedType,
  getMergedRenderableSchema,
  getMergedShaderCode,
} from './renderable/merged.js';

const getNextId = idFactory(0, 0x7fffffff);

export const getNextMaterialId = idFactory(0, 0x7fffffff);

export interface GraphicsRenderObject<T extends RenderObjectType = RenderObjectType> {
  readonly id: number;
  readonly type: T;
  readonly values: RenderObjectValues<T>;
  readonly state: RenderableState;
  readonly materialId: number;
}

export type RenderObjectType =
  | 'mesh'
  | 'points'
  | 'spheres'
  | 'cylinders'
  | 'text'
  | 'lines'
  | 'direct-volume'
  | 'image'
  | 'texture-mesh';

export type RenderObjectValues<T extends RenderObjectType> = T extends 'mesh'
  ? MeshValues
  : T extends 'points'
    ? PointsValues
    : T extends 'spheres'
      ? SpheresValues
      : T extends 'cylinders'
        ? CylindersValues
        : T extends 'text'
          ? TextValues
          : T extends 'lines'
            ? LinesValues
            : T extends 'direct-volume'
              ? DirectVolumeValues
              : T extends 'image'
                ? ImageValues
                : T extends 'texture-mesh'
                  ? TextureMeshValues
                  : never;

//

export function createRenderObject<T extends RenderObjectType>(
  type: T,
  values: RenderObjectValues<T>,
  state: RenderableState,
  materialId: number,
): GraphicsRenderObject<T> {
  return { id: getNextId(), type, values, state, materialId } as GraphicsRenderObject<T>;
}

//

/**
 * A render object that shares a single render item between the values of
 * multiple member render objects to reduce draw calls. The members are
 * not meant to be added to a scene themselves.
 */
export interface MergedGraphicsRenderObject<T extends RenderObjectType = RenderObjectType>
  extends GraphicsRenderObject<T> {
  readonly members: ReadonlyArray<GraphicsRenderObject<T>>;
  readonly merged: Merged;
}

export function isMergedRenderObject<T extends RenderObjectType>(
  o: GraphicsRenderObject<T>,
): o is MergedGraphicsRenderObject<T> {
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
  if (
    !canMergeValues(
      type,
      members.map((m) => m.values as MergeableValues),
    )
  )
    return false;
  return true;
}

/**
 * Creates a merged render object from member render objects, assumes
 * `canMergeRenderObjects(members)`. Shares the state of the first member.
 */
export function createMergedRenderObject<T extends RenderObjectType>(
  members: readonly GraphicsRenderObject<T>[],
): MergedGraphicsRenderObject<T> {
  if (!isMergeableType(members[0].type)) throw new Error(`unsupported merged render object type '${members[0].type}'`);
  const merged = createMergedValues(
    members[0].type,
    members.map((m) => m.values as MergeableValues),
  );
  return {
    id: getNextId(),
    type: members[0].type,
    values: merged.values as unknown as RenderObjectValues<T>,
    state: members[0].state,
    materialId: members[0].materialId,
    members,
    merged,
  };
}

function isTexture(x: unknown): x is Texture {
  return typeof x === 'object' && x !== null && typeof (x as Texture).attachFramebuffer === 'function';
}

/** The value cells of `o` (and of its members) holding textures not owned by `graphics` */
function getForeignTextureCells(o: GraphicsRenderObject, graphics: GraphicsContext): ValueCell<Texture>[] {
  const cells: ValueCell<Texture>[] = [];
  const add = (values: object) => {
    for (const k in values) {
      const cell = (values as Record<string, ValueCell<unknown> | undefined>)[k];
      const v = cell?.ref?.value;
      if (isTexture(v) && !isNullTexture(v) && !graphics.resources.ownsTexture(v)) {
        if (!cells.includes(cell as ValueCell<Texture>)) cells.push(cell as ValueCell<Texture>);
      }
    }
  };
  if (isMergedRenderObject(o)) for (const m of o.members) add(m.values);
  add(o.values);
  return cells;
}

/** Whether values of `o` hold textures of another graphics context than `graphics` */
export function hasForeignTextures(o: GraphicsRenderObject, graphics: GraphicsContext): boolean {
  return getForeignTextureCells(o, graphics).length > 0;
}

/**
 * Replaces textures of another context in the values of `o` by textures of `graphics` with the same contents, where
 * they were loaded from data (see `transferTexture`); false if some were rendered or computed (e.g. GPU surfaces),
 * which their producers create again
 */
export function transferTextures(o: GraphicsRenderObject, graphics: GraphicsContext): boolean {
  let complete = true;
  for (const cell of getForeignTextureCells(o, graphics)) {
    const texture = transferTexture(cell.ref.value, graphics);
    if (texture) ValueCell.update(cell, texture);
    else complete = false;
  }
  return complete;
}

/** What the programs of a renderable are made of: the kind of its render object, and whether it is merged */
export interface RenderableKind {
  readonly type: RenderObjectType;
  /** a merged render object (`SegmentSchema` added) */
  readonly merged: boolean;
}

const RenderablePrograms: { [T in RenderObjectType]: { shaderCode: ShaderCode; schema: RenderableSchema } } = {
  mesh: { shaderCode: MeshShaderCode, schema: MeshRenderableSchema },
  points: { shaderCode: PointsShaderCode, schema: PointsRenderableSchema },
  spheres: { shaderCode: SpheresShaderCode, schema: SpheresRenderableSchema },
  cylinders: { shaderCode: CylindersShaderCode, schema: CylindersRenderableSchema },
  text: { shaderCode: TextShaderCode, schema: TextRenderableSchema },
  lines: { shaderCode: LinesShaderCode, schema: LinesRenderableSchema },
  'direct-volume': { shaderCode: DirectVolumeShaderCode, schema: DirectVolumeRenderableSchema },
  image: { shaderCode: ImageShaderCode, schema: ImageRenderableSchema },
  'texture-mesh': { shaderCode: MeshShaderCode, schema: TextureMeshRenderableSchema },
};

/**
 * The shader code and schema the renderables of `kind` create their programs from (the same objects), `undefined`
 * for a kind that can't be merged
 */
export function getRenderableProgramSource(kind: RenderableKind) {
  if (!kind.merged) return RenderablePrograms[kind.type];
  if (!isMergeableType(kind.type)) return;
  const type = kind.type as MergedType;
  return { shaderCode: getMergedShaderCode(type), schema: getMergedRenderableSchema(type) };
}

const RenderableKinds = new WeakMap<Renderable<any>, RenderableKind>();

/** The kind `renderable` was created as, see `createRenderable` */
export function getRenderableKind(renderable: Renderable<any>): RenderableKind | undefined {
  return RenderableKinds.get(renderable);
}

export function createRenderable<T extends RenderObjectType>(
  ctx: GraphicsContext,
  o: GraphicsRenderObject<T>,
  transparency: Transparency,
  globals: GlobalDefines,
): Renderable<any> {
  const renderable = _createRenderable(ctx, o, transparency, globals);
  RenderableKinds.set(renderable, {
    type: isMergedRenderObject(o) ? o.merged.type : o.type,
    merged: isMergedRenderObject(o),
  });
  return renderable;
}

function _createRenderable<T extends RenderObjectType>(
  ctx: GraphicsContext,
  o: GraphicsRenderObject<T>,
  transparency: Transparency,
  globals: GlobalDefines,
): Renderable<any> {
  if (isMergedRenderObject(o)) {
    return MergedRenderable(
      ctx,
      o.id,
      o.merged,
      o.members.map((m) => m.values as MergeableValues),
      o.state,
      o.materialId,
      transparency,
      globals,
    );
  }
  switch (o.type) {
    case 'mesh':
      return MeshRenderable(ctx, o.id, o.values as MeshValues, o.state, o.materialId, transparency, globals);
    case 'points':
      return PointsRenderable(ctx, o.id, o.values as PointsValues, o.state, o.materialId, transparency, globals);
    case 'spheres':
      return SpheresRenderable(ctx, o.id, o.values as SpheresValues, o.state, o.materialId, transparency, globals);
    case 'cylinders':
      return CylindersRenderable(ctx, o.id, o.values as CylindersValues, o.state, o.materialId, transparency, globals);
    case 'text':
      return TextRenderable(ctx, o.id, o.values as TextValues, o.state, o.materialId, transparency, globals);
    case 'lines':
      return LinesRenderable(ctx, o.id, o.values as LinesValues, o.state, o.materialId, transparency, globals);
    case 'direct-volume':
      return DirectVolumeRenderable(
        ctx,
        o.id,
        o.values as DirectVolumeValues,
        o.state,
        o.materialId,
        transparency,
        globals,
      );
    case 'image':
      return ImageRenderable(ctx, o.id, o.values as ImageValues, o.state, o.materialId, transparency, globals);
    case 'texture-mesh':
      return TextureMeshRenderable(
        ctx,
        o.id,
        o.values as TextureMeshValues,
        o.state,
        o.materialId,
        transparency,
        globals,
      );
  }
  assertUnreachable(o.type);
}
