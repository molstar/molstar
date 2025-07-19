/**
 * Copyright (c) 2021 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Michal Malý <malym@ibt.cas.cz>
 */

const fs = require('fs');
const path = require('path');
const argparse = require('argparse');

function removeDir(dirPath) {
    for (const ent of fs.readdirSync(dirPath)) {
        const entryPath = path.join(dirPath, ent);
        remove(entryPath);
    }

    fs.rmdirSync(dirPath);
}

function remove(entryPath) {
    const st = fs.statSync(entryPath);
    if (st.isDirectory())
        removeDir(entryPath);
    else
        fs.unlinkSync(entryPath);
}

const argParser = new argparse.ArgumentParser({
    add_help: true,
    description: 'Clean Script'
});
argParser.add_argument('--build', { required: false, action: 'store_true' });
argParser.add_argument('--lib', { required: false, action: 'store_true' });
argParser.add_argument('--all', { required: false, action: 'store_true' });
const args = argParser.parse_args();

const toClean = [];

if (args.build || args.all) {
    toClean.push(path.resolve(__dirname, '../build'));
    toClean.push(path.resolve(__dirname, '../deploy/data'));
}
if (args.lib || args.all) {
    toClean.push(
        path.resolve(__dirname, '../lib'),
        path.resolve(__dirname, '../tsconfig.tsbuildinfo'),
    );
}

console.log('\n###', 'cleaning', toClean.join(', '));

toClean.forEach(ph => {
    if (fs.existsSync(ph)) {
        try {
            remove(ph);
        } catch (err) {
            console.warn(`Cleanup failed: ${err}`);
        }
    }
});
