/**
 * Copyright (c) 2017-2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import * as argparse from 'argparse'
import * as fs from 'fs'

import { validate } from './util/validate'

function runValidateSchema (path: string) {
    const str = fs.readFileSync(path, 'utf8')
    const result = validate(JSON.parse(str))
    console.log(result === true ? 'valid json schema' : `invalid json schema: "${result}"`)
}

const parser = new argparse.ArgumentParser({
    addHelp: true,
    description: 'Validate json schema'
});
parser.addArgument([ 'path' ], {
    help: 'path to json schema'
});
const args = parser.parseArgs();

if (args.path) {
    runValidateSchema(args.path)
}
