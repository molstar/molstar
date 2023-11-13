/**
 * Copyright (c) 2023 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Adam Midlik <midlik@gmail.com>
 *
 * Command-line application for rendering images from MolViewSpec files
 * Build: npm install --no-save gl jpeg-js pngjs  // these packages are not listed in dependencies for performance reasons
 *        npm run build
 * Run:   node lib/commonjs/examples/mvs-render -i examples/1cbs.mvsj -o ../outputs/1cbs.png
 */

import { ArgumentParser } from 'argparse';
import fs from 'fs';
import gl from 'gl';
import jpegjs from 'jpeg-js';
import path from 'path';
import pngjs from 'pngjs';

import { MolViewSpec } from '../../extensions/mvs/behavior';
import { loadMVS } from '../../extensions/mvs/load';
import { MVSData } from '../../extensions/mvs/mvs-data';
import { HeadlessPluginContext } from '../../mol-plugin/headless-plugin-context';
import { DefaultPluginSpec, PluginSpec } from '../../mol-plugin/spec';
import { ExternalModules } from '../../mol-plugin/util/headless-screenshot';
import { setFSModule } from '../../mol-util/data-source';
import { treeToString } from '../../extensions/mvs/tree/generic/tree-utils';


setFSModule(fs);

const DEFAULT_SIZE = '800x800';

interface Args {
    input: string[],
    output: string[],
    size: { width: number, height: number },
}

function parseArguments(): Args {
    const parser = new ArgumentParser({ description: 'Command-line application for rendering images from MolViewSpec files' });
    parser.add_argument('-i', '--input', { required: true, nargs: '+', help: 'Input file(s) in .mvsj format' });
    parser.add_argument('-o', '--output', { required: true, nargs: '+', help: 'File path(s) for output files (one output path for each input file). Output format is inferred from the file extension (.png or .jpg)' });
    parser.add_argument('-s', '--size', { help: `Output image resolution, {width}x{height}. Default: ${DEFAULT_SIZE}.`, default: DEFAULT_SIZE });
    const args = parser.parse_args();
    try {
        const parts = args.size.split('x');
        if (parts.length !== 2) throw new Error('Must contain two x-separated parts');
        args.size = { width: parseIntStrict(parts[0]), height: parseIntStrict(parts[1]) };
    } catch {
        parser.error(`argument: --size: invalid image size string: '${args.size}' (must be two x-separated integers (width and height), e.g. '400x300')`);
    }
    if (args.input.length !== args.output.length) {
        parser.error(`argument: --output: must specify the same number of input and output file paths (specified ${args.input.length} input paths, but ${args.output.length} output paths)`);
    }
    return { ...args };
}

async function main() {
    const args = parseArguments();
    console.log('input:', args.input);
    console.log('output:', args.output);

    // Create a headless plugin
    const externalModules: ExternalModules = { gl, pngjs, 'jpeg-js': jpegjs };
    const spec = DefaultPluginSpec();
    spec.behaviors.push(PluginSpec.Behavior(MolViewSpec));
    const plugin = new HeadlessPluginContext(externalModules, spec, args.size);
    await plugin.init();
    console.log('init done')

    for (let i = 0; i < args.input.length; i++) {
        const input = args.input[i];
        const output = args.output[i];
        const data = fs.readFileSync(input, { encoding: 'utf8' });
        const mvsData = MVSData.fromMVSJ(data);
        console.log('read done', treeToString(mvsData.root))

        await loadMVS(plugin, mvsData, { sanityChecks: true, deletePrevious: false });
        console.log('load done', input)
        fs.mkdirSync(path.dirname(output), { recursive: true });
        await plugin.saveStateSnapshot(output + '.molj')
        console.log('save molj done', output)
        await plugin.saveImage(output, undefined, undefined);
        console.log('save done', output)
    }
    await plugin.clear();
    console.log('clear done')
    plugin.dispose();
}

/** Parse integer, fail early. */
function parseIntStrict(str: string): number {
    if (str === '') throw new Error('Is empty string');
    const result = Number(str);
    if (isNaN(result)) throw new Error('Is NaN');
    if (Math.floor(result) !== result) throw new Error('Is not integer');
    return result;
}

main();
