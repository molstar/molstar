/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Taken/adapted from DensityServer (https://github.com/dsehnal/DensityServer)
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import pack from './pack/main'
import VERSION from './pack/version'

interface Config {
    input: { name: string, filename: string }[],
    format: 'ccp4' | 'dsn6',
    isPeriodic: boolean,
    outputFilename: string,
    blockSizeInMB: number
}

let config: Config = {
    input: [],
    format: 'ccp4',
    isPeriodic: false,
    outputFilename: '',
    blockSizeInMB: 96
};

function getFormat(format: string): Config['format'] {
    switch (format.toLowerCase()) {
        case 'ccp4': return 'ccp4'
        case 'dsn6': return 'dsn6'
    }
    throw new Error(`unsupported format '${format}'`)
}

function printHelp() {
    let help = [
        `VolumeServer Packer ${VERSION}, (c) 2016 - now, David Sehnal`,
        ``,
        `The input data must be CCP4/MAP mode 2 (32-bit floats) files.`,
        ``,
        `Usage: `,
        ``,
        `  node pack -v`,
        `    Print version.`,
        ``,
        `  node pack -xray main.ccp4 diff.ccp4 output.mdb [-blockSize 96]`,
        `    Pack main and diff density into a single block file.`,
        `    Optionally specify maximum block size.`,
        ``,
        `  node pack -em density.map output.mdb [-blockSize 96]`,
        `    Pack single density into a block file.`,
        `    Optionally specify maximum block size.`
    ];
    console.log(help.join('\n'));
}

function parseInput() {
    let input = false;

    if (process.argv.length <= 2) {
        printHelp();
        process.exit();
        return false;
    }

    for (let i = 2; i < process.argv.length; i++) {
        switch (process.argv[i].toLowerCase()) {
            case '-blocksize':
                config.blockSizeInMB = +process.argv[++i];
                break;
            case '-format':
                config.format = getFormat(process.argv[++i]);
                break;
            case '-xray':
                input = true;
                config.input = [
                    { name: '2Fo-Fc', filename: process.argv[++i] },
                    { name: 'Fo-Fc', filename: process.argv[++i] }
                ];
                config.isPeriodic = true;
                config.outputFilename = process.argv[++i];
                break;
            case '-em':
                input = true;
                config.input = [
                    { name: 'em', filename: process.argv[++i] }
                ];
                config.outputFilename = process.argv[++i];
                break;
            case '-v':
                console.log(VERSION);
                process.exit();
                return false;
            default:
                printHelp();
                process.exit();
                return false;
        }
    }
    return input;
}

if (parseInput()) {
    pack(config.input, config.blockSizeInMB, config.isPeriodic, config.outputFilename, config.format);
}