/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import fetch from 'node-fetch';
import * as fs from 'fs'
import * as path from 'path'
import * as argparse from 'argparse'
import { makeDir } from 'mol-util/make-dir';
import { now } from 'mol-task';
import { PerformanceMonitor } from 'mol-util/performance-monitor';

const cmdParser = new argparse.ArgumentParser({
    addHelp: true,
    description: 'Download JSON data from PDBe API'
});

cmdParser.addArgument(['--in'], { help: 'Input folder', required: true });
cmdParser.addArgument(['--out'], { help: 'Output folder', required: true });

interface CmdArgs {
    in: string,
    out: string
}

const cmdArgs = cmdParser.parseArgs() as CmdArgs;

function getPDBid(name: string) {
    let idx = name.indexOf('_');
    if (idx < 0) idx = name.indexOf('.');
    return name.substr(0, idx).toLowerCase();
}

function findEntries() {
    const files = fs.readdirSync(cmdArgs.in);
    const cifTest = /\.cif$/;
    const groups = new Map<string, string[]>();
    const keys: string[] = [];

    for (const f of files) {
        if (!cifTest.test(f)) continue;
        const id = getPDBid(f);
        const groupId = `${id[1]}${id[2]}`;

        if (groups.has(groupId)) groups.get(groupId)!.push(id);
        else {
            keys.push(groupId);
            groups.set(groupId, [id]);
        }
    }

    const ret: { key: string, entries: string[] }[] = [];
    for (const key of keys) {
        ret.push({ key, entries: groups.get(key)! })
    }

    return ret;
}

async function process() {
    const entries = findEntries();
    makeDir(cmdArgs.out);

    const started = now();
    let prog = 0;
    for (const e of entries) {
        const ts = now();
        console.log(`${prog}/${entries.length} ${e.entries.length} entries.`)
        const body = e.entries.join(',');
        const query = await fetch(`https://www.ebi.ac.uk/pdbe/api/validation/residuewise_outlier_summary/entry`, { method: 'POST', body });
        const data = await query.text();
        fs.writeFileSync(path.join(cmdArgs.out, e.key + '.json'), data);
        const time = now() - started;
        console.log(`${++prog}/${entries.length} in ${PerformanceMonitor.format(time)} (last ${PerformanceMonitor.format(now() - ts)}, avg ${PerformanceMonitor.format(time / prog)})`);
    }
}

process();