#!/usr/bin/env node
/**
 * Copyright (c) 2023 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Adam Midlik <midlik@gmail.com>
 *
 * Command-line application for printing MolViewSpec tree schema
 * Build: npm run build
 * Run:   node lib/commonjs/cli/mvs/mvs-print-schema
 *        node lib/commonjs/cli/mvs/mvs-print-schema --markdown
 */

import { ArgumentParser } from 'argparse';
import { treeSchemaToMarkdown, treeSchemaToString } from '../../extensions/mvs/tree/generic/tree-schema';
import { MVSTreeSchema } from '../../extensions/mvs/tree/mvs/mvs-tree';


/** Command line argument values for `main` */
interface Args {
    markdown: boolean,
}

/** Return parsed command line arguments for `main` */
function parseArguments(): Args {
    const parser = new ArgumentParser({ description: 'Command-line application for printing MolViewSpec tree schema.' });
    parser.add_argument('-m', '--markdown', { action: 'store_true', help: 'Print the schema as markdown instead of plain text.' });
    const args = parser.parse_args();
    return { ...args };
}

/** Main workflow for printing MolViewSpec tree schema. */
function main(args: Args) {
    if (args.markdown) {
        console.log(treeSchemaToMarkdown(MVSTreeSchema));
    } else {
        console.log(treeSchemaToString(MVSTreeSchema));
    }
}

main(parseArguments());
