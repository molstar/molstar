/**
 * Copyright (c) 2023 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Adam Midlik <midlik@gmail.com>
 *
 * Command-line application for rendering images from MolViewSpec files
 * Build: npm install --no-save gl jpeg-js pngjs  // these packages are not listed in Mol* dependencies for performance reasons
 *        npm run build
 * Run:   node lib/commonjs/examples/mvs/mvs-render -i examples/mvs/1cbs.mvsj -o ../outputs/1cbs.png --size 800x600 --molj
 */

import { ArgumentParser } from 'argparse';
import fs from 'fs';
import gl from 'gl';
import jpegjs from 'jpeg-js';
import path from 'path';
import pngjs from 'pngjs';

import { Canvas3DParams } from '../../mol-canvas3d/canvas3d';
import { PluginContext } from '../../mol-plugin/context';
import { HeadlessPluginContext } from '../../mol-plugin/headless-plugin-context';
import { DefaultPluginSpec, PluginSpec } from '../../mol-plugin/spec';
import { ExternalModules, defaultCanvas3DParams } from '../../mol-plugin/util/headless-screenshot';
import { setFSModule } from '../../mol-util/data-source';
import { onelinerJsonString } from '../../mol-util/json';
import { ParamDefinition as PD } from '../../mol-util/param-definition';

// MolViewSpec must be imported after HeadlessPluginContext
import { MolViewSpec } from '../../extensions/mvs/behavior';
import { loadMVS } from '../../extensions/mvs/load';
import { MVSData } from '../../extensions/mvs/mvs-data';
import { dfs } from '../../extensions/mvs/tree/generic/tree-utils';


setFSModule(fs);

const DEFAULT_SIZE = '800x800';

/** Command line argument values for `main` */
interface Args {
    input: string[],
    output: string[],
    size: { width: number, height: number },
    molj: boolean,
}

/** Return parsed command line arguments for `main` */
function parseArguments(): Args {
    const parser = new ArgumentParser({ description: 'Command-line application for rendering images from MolViewSpec files' });
    parser.add_argument('-i', '--input', { required: true, nargs: '+', help: 'Input file(s) in .mvsj format' });
    parser.add_argument('-o', '--output', { required: true, nargs: '+', help: 'File path(s) for output files (one output path for each input file). Output format is inferred from the file extension (.png or .jpg)' });
    parser.add_argument('-s', '--size', { help: `Output image resolution, {width}x{height}. Default: ${DEFAULT_SIZE}.`, default: DEFAULT_SIZE });
    parser.add_argument('-m', '--molj', { action: 'store_true', help: `Save Mol* state (.molj) in addition to rendered images (use the same output file paths but with .molj extension)` });
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

        const data = fs.readFileSync(input, { encoding: 'utf8' });
        const mvsData = MVSData.fromMVSJ(data);
        removeLabelNodes(mvsData);

        await loadMVS(plugin, mvsData, { sanityChecks: true, replaceExisting: true });
        fs.mkdirSync(path.dirname(output), { recursive: true });
        if (args.molj) {
            await plugin.saveStateSnapshot(withExtension(output, '.molj'));
        }
        await plugin.saveImage(output);
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

/** Remove any label* nodes from the MVS tree (in-place). Print warning if removed at least one node. */
function removeLabelNodes(mvsData: MVSData): void {
    let removedSomething = false;
    dfs(mvsData.root, node => {
        const nChildren = node.children?.length ?? 0;
        node.children = node.children?.filter(c => !c.kind.startsWith('label'));
        if ((node.children?.length ?? 0) !== nChildren) {
            removedSomething = true;
        }
    });
    if (removedSomething) {
        // trying to render labels would fail because `document` is not available in Nodejs
        console.error('Rendering labels is not yet supported in mvs-render. Skipping any label* nodes.');
    }
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
