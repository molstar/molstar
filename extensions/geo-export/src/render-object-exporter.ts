/**
 * Copyright (c) 2021 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Sukolsak Sakshuwong <sukolsak@stanford.edu>
 */

import type { GraphicsRenderObject } from '@molstar/graphics/gl/render-object';
import type { GraphicsContext } from '@molstar/graphics/gl/context';
import { RuntimeContext } from '@molstar/core/task';

export type RenderObjectExportData = {
  [k: string]: string | Uint8Array | ArrayBuffer | undefined;
};

export interface RenderObjectExporter<D extends RenderObjectExportData> {
  readonly fileExtension: string;
  add(renderObject: GraphicsRenderObject, graphics: GraphicsContext, ctx: RuntimeContext): Promise<void> | undefined;
  getData(ctx: RuntimeContext): Promise<D>;
  getBlob(ctx: RuntimeContext): Promise<Blob>;
}
