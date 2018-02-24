/**
 * Copyright (c) 2017-2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import * as argparse from 'argparse'
import * as fs from 'fs'

import { generate } from './util/generate'

function generateSchema (name: string, path: string) {
    const str = fs.readFileSync(path, 'utf8')
    return generate(name, JSON.parse(str))
}

const parser = new argparse.ArgumentParser({
    addHelp: true,
    description: 'Argparse example'
});
parser.addArgument([ 'name' ], {
    help: 'schema name'
});
parser.addArgument([ 'path' ], {
    help: 'json schema file path'
});
parser.addArgument([ '--out', '-o' ], {
    help: 'generated typescript output path'
});
const args = parser.parseArgs();

if (args.name && args.path) {
    const schema = generateSchema(args.name, args.path)
    if (args.out) {
        fs.writeFileSync(args.out, schema)
    } else {
        console.log(schema)
    }
}
