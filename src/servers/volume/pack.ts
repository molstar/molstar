#!/usr/bin/env node
/**
 * Copyright (c) 2018-2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Taken/adapted from DensityServer (https://github.com/dsehnal/DensityServer)
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import * as argparse from 'argparse';
import { pack } from './pack/main';
import { VERSION } from './pack/version';

type FileFormat = 'ccp4' | 'dsn6'

interface Config {
    input: { name: string, filename: string }[],
    format: FileFormat,
    isPeriodic: boolean,
    outputFilename: string,
    blockSizeInMB: number
}

function getConfig(args: Args) {
    const config: Partial<Config> = {
        blockSizeInMB: args.blockSizeInMB,
        format: args.format,
        outputFilename: args.output
    };
    switch (args.mode) {
        case 'em':
            config.input = [
                { name: 'em', filename: args.inputEm }
            ];
            config.isPeriodic = false;
            break;
        case 'xray':
            config.input = [
                { name: '2Fo-Fc', filename: args.input2fofc },
                { name: 'Fo-Fc', filename: args.inputFofc }
            ];
            config.isPeriodic = true;
            break;
    }
    return config as Config;
}

interface GeneralArgs {
    blockSizeInMB: number
    format: FileFormat
    output: string
}
interface XrayArgs extends GeneralArgs {
    mode: 'xray'
    input2fofc: string
    inputFofc: string
}
interface EmArgs extends GeneralArgs {
    mode: 'em'
    inputEm: string
}
type Args = XrayArgs | EmArgs

const parser = new argparse.ArgumentParser({
    add_help: true,
    description: `VolumeServer Packer ${VERSION}, (c) 2018-2019, Mol* contributors`
});

const subparsers = parser.add_subparsers({
    title: 'Packing modes',
    dest: 'mode'
});

function addGeneralArgs(parser: argparse.ArgumentParser) {
    parser.add_argument('output', { help: `Output path.` });
    parser.add_argument('--blockSizeInMB', { default: 96, help: `Maximum block size.`, metavar: 'SIZE' });
    parser.add_argument('--format', { default: 'ccp4', help: `Input file format.` });
}

const xrayParser = subparsers.add_parser('xray', { add_help: true });
xrayParser.add_argument('input2fofc', { help: `Path to 2fofc file.`, metavar: '2FOFC' });
xrayParser.add_argument('inputFofc', { help: `Path to fofc file.`, metavar: 'FOFC' });
addGeneralArgs(xrayParser);

const emParser = subparsers.add_parser('em', { add_help: true });
emParser.add_argument('inputEm', { help: `Path to EM density file.`, metavar: 'EM' });
addGeneralArgs(emParser);

const args: Args = parser.parse_args();
const config = getConfig(args);

pack(config.input, config.blockSizeInMB, config.isPeriodic, config.outputFilename, config.format);
