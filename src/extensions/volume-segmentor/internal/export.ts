/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Tadej Satler <tadej.satler@gmail.com>
 */

import { Volume } from '../../../mol-model/volume';
import { RuntimeContext } from '../../../mol-task';
import { CCP4Writer } from '../../../mol-io/writer/ccp4/ccp4';
import { download } from '../../../mol-util/download';
import { zip } from '../../../mol-util/zip/zip';
import { BodyInfo, BodyMaskParams, BodyMaskResult, LabelStore } from '../types';
import { computeBodyMask, scatterToFullBox } from './mask-compute';

/** Bodies with voxels, largest first (ties broken by id) — the order masks are exported in. */
export function orderBodiesBySize(bodies: readonly BodyInfo[]): BodyInfo[] {
    return bodies
        .filter(b => b.voxelCount > 0)
        .sort((a, b) => b.voxelCount - a.voxelCount || a.id - b.id);
}

/** File-name stem of a volume label (extension stripped), falling back to `volume`. */
export function maskBaseName(label: string | undefined): string {
    const base = (label ?? '').replace(/\.gz$/i, '').replace(/\.[^.]+$/, '').trim();
    return base || 'volume';
}

/** `{base}_body001_mask.mrc`, `rank` is 1-based. */
export function bodyMaskFileName(base: string, rank: number): string {
    return `${base}_body${String(rank).padStart(3, '0')}_mask.mrc`;
}

/** Effective mask params of one body: per-body overrides on top of the defaults. */
export function resolveBodyMaskParams(body: BodyInfo, defaults: BodyMaskParams): BodyMaskParams {
    return {
        extend: body.extend ?? defaults.extend,
        softEdge: body.softEdge ?? defaults.softEdge,
        pruneBelowThreshold: defaults.pruneBelowThreshold,
    };
}

/** MRC bytes of a body mask on the full source grid (same box and origin as the source). */
export function writeBodyMaskMrc(volume: Volume, result: BodyMaskResult): ArrayBuffer {
    const full = new Float32Array(volume.grid.cells.data.length);
    scatterToFullBox(result, volume.grid.cells.space, full);
    return CCP4Writer.writeMrc(volume.grid, full);
}

/** Downloads the volume's current voxel data as a float32 MRC, in-place edits included. */
export function downloadVolumeMrc(volume: Volume, baseName: string) {
    const raw = volume.grid.cells.data as unknown as ArrayLike<number>;
    const data = raw instanceof Float32Array ? raw : Float32Array.from(raw);
    const buffer = CCP4Writer.writeMrc(volume.grid, data);
    download(new Blob([buffer], { type: 'application/octet-stream' }), `${baseName}_volume.mrc`);
}

export interface ExportBodyMasksOptions {
    baseName: string
    /** Bundle all masks into one zip instead of downloading them one by one. */
    bundleZip: boolean
    /** Deflate the zip (slow for large maps); ignored unless `bundleZip`. */
    compress: boolean
}

/**
 * Computes and downloads one MRC mask per body, largest body first. Returns the file names
 * in export order.
 */
export async function exportBodyMasks(volume: Volume, store: LabelStore, defaults: BodyMaskParams, thresholdAbs: number, options: ExportBodyMasksOptions, ctx: RuntimeContext): Promise<string[]> {
    const ordered = orderBodiesBySize(store.bodies);
    const files: { [name: string]: Uint8Array<ArrayBuffer> } = {};
    const names: string[] = [];

    for (let i = 0; i < ordered.length; i++) {
        const body = ordered[i];
        await ctx.update({ message: `Computing mask for ${body.name}…`, current: i, max: ordered.length });
        const result = await computeBodyMask(volume, store.labels, body.id, resolveBodyMaskParams(body, defaults), thresholdAbs, ctx);
        if (!result || result.voxelCount === 0) continue;

        const name = bodyMaskFileName(options.baseName, names.length + 1);
        names.push(name);
        const buffer = writeBodyMaskMrc(volume, result);
        if (options.bundleZip) {
            files[name] = new Uint8Array(buffer);
        } else {
            download(new Blob([buffer], { type: 'application/octet-stream' }), name);
        }
    }

    if (options.bundleZip && names.length > 0) {
        await ctx.update({ message: 'Bundling masks…' });
        const zipped = await zip(ctx, files, !options.compress);
        download(new Blob([zipped], { type: 'application/zip' }), `${options.baseName}_body_masks.zip`);
    }
    return names;
}
