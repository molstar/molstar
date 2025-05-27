#!/usr/bin/env node
/**
 * Copyright (c) 2023-2024 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Adam Midlik <midlik@gmail.com>
 *
 * Command-line application for rendering images from MolViewSpec files
 * From Molstar NPM package:
 *     npm install molstar canvas gl jpeg-js pngjs
 *     npx mvs-render -i examples/mvs/1cbs.mvsj -o ../outputs/1cbs.png --size 800x600 --molj
 * From Molstar source code:
 *     npm install
 *     npm install --no-save canvas gl jpeg-js pngjs  // these packages are not listed in Mol* dependencies for performance reasons
 *     npm run build
 *     node lib/commonjs/cli/mvs/mvs-render -i examples/mvs/1cbs.mvsj -o ../outputs/1cbs.png --size 800x600 --molj
 */

import { ArgumentParser } from 'argparse';
import fs from 'fs';
import gl from 'gl';
import jpegjs from 'jpeg-js';
import path from 'path';
import pngjs from 'pngjs';

import { Canvas3DParams } from '../../mol-canvas3d/canvas3d';
import { setCanvasModule } from '../../mol-geo/geometry/text/font-atlas';
import { PluginContext } from '../../mol-plugin/context';
import { HeadlessPluginContext } from '../../mol-plugin/headless-plugin-context';
import { DefaultPluginSpec, PluginSpec } from '../../mol-plugin/spec';
import { ExternalModules, defaultCanvas3DParams } from '../../mol-plugin/util/headless-screenshot';
import { Task } from '../../mol-task';
import { setFSModule } from '../../mol-util/data-source';
import { onelinerJsonString } from '../../mol-util/json';
import { ParamDefinition as PD } from '../../mol-util/param-definition';

// MolViewSpec must be imported after HeadlessPluginContext
import { Mp4Export } from '../../extensions/mp4-export';
import { MolViewSpec } from '../../extensions/mvs/behavior';
import { loadMVSX } from '../../extensions/mvs/components/formats';
import { loadMVS } from '../../extensions/mvs/load';
import { MVSData } from '../../extensions/mvs/mvs-data';


setFSModule(fs);
setCanvasModule(require('canvas'));

const DEFAULT_SIZE = '800x800';

/** Command line argument values for `main` */
interface Args {
    input: string[],
    output: string[],
    size: { width: number, height: number },
    molj: boolean,
    no_extensions: boolean,
}

/** Return parsed command line arguments for `main` */
function parseArguments(): Args {
    const parser = new ArgumentParser({ description: 'Command-line application for rendering images from MolViewSpec files' });
    parser.add_argument('-i', '--input', { required: true, nargs: '+', help: 'Input file(s) in .mvsj or .mvsx format. File format is inferred from the file extension.' });
    parser.add_argument('-o', '--output', { required: true, nargs: '+', help: 'File path(s) for output files (one output path for each input file). Output format is inferred from the file extension (.png or .jpg)' });
    parser.add_argument('-s', '--size', { help: `Output image resolution, {width}x{height}. Default: ${DEFAULT_SIZE}.`, default: DEFAULT_SIZE });
    parser.add_argument('-m', '--molj', { action: 'store_true', help: `Save Mol* state (.molj) in addition to rendered images (use the same output file paths but with .molj extension)` });
    parser.add_argument('-n', '--no-extensions', { action: 'store_true', help: `Do not apply builtin MVS-loading extensions (not a part of standard MVS specification)` });
    const args = parser.parse_args();
    try {
        const parts = args.size.split('x');
        if (parts.length !== 2) throw new Error('Must contain two x-separated parts');
        args.size = { width: parseIntStrict(parts[0]), height: parseIntStrict(parts[1]) };
    } catch {
        parser.error(`argument: --size: invalid image size string: '${args.size}' (must be two x-separated integers (width and height), e.g. '400x300')`);
    }
    if (args.input.length !== args.output.length) {
        parser.error(`argument: --output: must specify the same number of input and output file paths (specified ${args.input.length} input path${args.input.length !== 1 ? 's' : ''} but ${args.output.length} output path${args.output.length !== 1 ? 's' : ''})`);
    }
    return { ...args };
}

/** Main workflow for rendering images from MolViewSpec files */
async function main(args: Args): Promise<void> {
    const plugin = await createHeadlessPlugin(args);

    for (let i = 0; i < args.input.length; i++) {
        const input = args.input[i];
        const output = args.output[i];
        console.log(`Processing ${input} -> ${output}`);

        let mvsData: MVSData;
        let sourceUrl: string | undefined;
        if (input.toLowerCase().endsWith('.mvsj')) {
            const data = fs.readFileSync(input, { encoding: 'utf8' });
            mvsData = MVSData.fromMVSJ(data);
            sourceUrl = `file://${path.resolve(input)}`;
        } else if (input.toLowerCase().endsWith('.mvsx')) {
            const data = fs.readFileSync(input);
            const mvsx = await plugin.runTask(Task.create('Load MVSX', async ctx => loadMVSX(plugin, ctx, data)));
            mvsData = mvsx.mvsData;
            sourceUrl = mvsx.sourceUrl;
        } else {
            throw new Error(`Input file name must end with .mvsj or .mvsx: ${input}`);
        }
        await loadMVS(plugin, mvsData, { sanityChecks: true, sourceUrl: sourceUrl, extensions: args.no_extensions ? [] : undefined });

        fs.mkdirSync(path.dirname(output), { recursive: true });
        if (args.molj) {
            await plugin.saveStateSnapshot(withExtension(output, '.molj'));
        }
        if (output.toLowerCase().endsWith('.mp4')) {
            await plugin.saveAnimation(output);
        } else {
            await plugin.saveImage(output);
        }
        checkState(plugin);
    }
    await plugin.clear();
    plugin.dispose();
}

/** Return a new and initiatized HeadlessPlugin */
async function createHeadlessPlugin(args: Pick<Args, 'size'>): Promise<HeadlessPluginContext> {
    const externalModules: ExternalModules = { gl, pngjs, 'jpeg-js': jpegjs };
    const spec = DefaultPluginSpec();
    spec.behaviors.push(PluginSpec.Behavior(MolViewSpec));
    spec.behaviors.push(PluginSpec.Behavior(Mp4Export));
    const headlessCanvasOptions = defaultCanvas3DParams();
    const canvasOptions = {
        ...PD.getDefaultValues(Canvas3DParams),
        cameraResetDurationMs: headlessCanvasOptions.cameraResetDurationMs,
        postprocessing: headlessCanvasOptions.postprocessing,
    };
    const plugin = new HeadlessPluginContext(externalModules, spec, args.size, { canvas: canvasOptions });
    try {
        await plugin.init();
    } catch (error) {
        plugin.dispose();
        throw error;
    }
    return plugin;
}

/** Parse integer, fail early. */
function parseIntStrict(str: string): number {
    if (str === '') throw new Error('Is empty string');
    const result = Number(str);
    if (isNaN(result)) throw new Error('Is NaN');
    if (Math.floor(result) !== result) throw new Error('Is not integer');
    return result;
}

/** Replace the file extension in `filename` by `extension`. If `filename` has no extension, add it. */
function withExtension(filename: string, extension: string): string {
    const oldExtension = path.extname(filename);
    return filename.slice(0, -oldExtension.length) + extension;
}

/** Check Mol* state, print and throw error if any cell is not OK. */
function checkState(plugin: PluginContext): void {
    const cells = Array.from(plugin.state.data.cells.values());
    const badCell = cells.find(cell => cell.status !== 'ok');
    if (badCell) {
        console.error(`Building Mol* state failed`);
        console.error(`    Transformer: ${badCell.transform.transformer.id}`);
        console.error(`    Params: ${onelinerJsonString(badCell.transform.params)}`);
        console.error(`    Error: ${badCell.errorText}`);
        console.error(``);
        throw new Error(`Building Mol* state failed: ${badCell.errorText}`);
    }
}

main(parseArguments());
