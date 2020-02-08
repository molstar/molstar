/**
 * Copyright (c) 2019-2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import * as express from 'express'
import * as compression from 'compression'
import * as cors from 'cors'
import * as bodyParser from 'body-parser'
import * as argparse from 'argparse'
import * as fs from 'fs'
import * as path from 'path'
import { makeDir } from '../../mol-util/make-dir'

interface Config {
    working_folder: string,
    port?: string | number,
    app_prefix: string,
    max_states: number
}

const cmdParser = new argparse.ArgumentParser({
    addHelp: true
});
cmdParser.addArgument(['--working-folder'], { help: 'Working forlder path.', required: true });
cmdParser.addArgument(['--port'], { help: 'Server port. Altenatively use ENV variable PORT.', type: 'int', required: false });
cmdParser.addArgument(['--app-prefix'], { help: 'Server app prefix.', defaultValue: '', required: false });
cmdParser.addArgument(['--max-states'], { help: 'Maxinum number of states that could be saved.', defaultValue: 40, type: 'int', required: false });

const Config = cmdParser.parseArgs() as Config;

if (!Config.port) Config.port = process.env.port || 1339;

const app = express();
app.use(compression(<any>{ level: 6, memLevel: 9, chunkSize: 16 * 16384, filter: () => true }));
app.use(cors({ methods: ['GET', 'PUT'] }));
app.use(bodyParser.json({ limit: '20mb' }));

type Index = { timestamp: number, id: string, name: string, description: string, isSticky?: boolean }[]

function createIndex() {
    const fn = path.join(Config.working_folder, 'index.json');
    if (fs.existsSync(fn)) return;
    if (!fs.existsSync(Config.working_folder)) makeDir(Config.working_folder);
    fs.writeFileSync(fn, '[]', 'utf-8');
}

function writeIndex(index: Index) {
    const fn = path.join(Config.working_folder, 'index.json');
    if (!fs.existsSync(Config.working_folder)) makeDir(Config.working_folder);
    fs.writeFileSync(fn, JSON.stringify(index, null, 2), 'utf-8');
}

function readIndex() {
    const fn = path.join(Config.working_folder, 'index.json');
    if (!fs.existsSync(fn)) return [];
    return JSON.parse(fs.readFileSync(fn, 'utf-8')) as Index;
}

function validateIndex(index: Index) {
    if (index.length > Config.max_states) {
        const deletes: Index = [], newIndex: Index = [];
        const toDelete = index.length - Config.max_states; 
        
        for (const e of index) {
            if (!e.isSticky && deletes.length < toDelete) {
                deletes.push(e);
            } else {
                newIndex.push(e);
            }
        }
        // index.slice(0, index.length - 30);
        for (const d of deletes) {
            try {
                fs.unlinkSync(path.join(Config.working_folder, d.id + '.json'))
            } catch { }
        }
        return newIndex;
    }
    return index;
}

function remove(id: string) {
    let index = readIndex();
    let i = 0;
    for (const e of index) {
        if (e.id !== id) {            
            i++;
            continue;
        }
        try {
            for (let j = i + 1; j < index.length; j++) {
                index[j - 1] = index[j];
            }
            index.pop();
            writeIndex(index);
        } catch { }
        try {
            fs.unlinkSync(path.join(Config.working_folder, e.id + '.json'))
        } catch { }
        return;
    }
}

function clear() {
    let index = readIndex();
    for (const e of index) {
        try {
            fs.unlinkSync(path.join(Config.working_folder, e.id + '.json'))
        } catch { }
    }
    writeIndex([]);
}

export function createv4() {
    let d = (+new Date());
    const uuid = 'xxxxxxxx-xxxx-4xxx-yxxx-xxxxxxxxxxxx'.replace(/[xy]/g, function(c) {
        const r = (d + Math.random()*16)%16 | 0;
        d = Math.floor(d/16);
        return (c==='x' ? r : (r&0x3|0x8)).toString(16);
    });
    return uuid.toLowerCase();
}

function mapPath(path: string) {
    if (!Config.app_prefix) return path;
    return `/${Config.app_prefix}/${path}`;
}

app.get(mapPath(`/get/:id`), (req, res) => {
    const id: string = req.params.id || '';
    console.log('Reading', id);
    if (id.length === 0 || id.indexOf('.') >= 0 || id.indexOf('/') >= 0 || id.indexOf('\\') >= 0) {
        res.status(404);
        res.end();
        return;
    }

    fs.readFile(path.join(Config.working_folder, id + '.json'), 'utf-8', (err, data) => {
        if (err) {
            res.status(404);
            res.end();
            return;
        }

        res.writeHead(200, {
            'Content-Type': 'application/json; charset=utf-8',
        });
        res.write(data);
        res.end();
    });
});

app.get(mapPath(`/clear`), (req, res) => {
    clear();
    res.status(200);
    res.end();
});

app.get(mapPath(`/remove/:id`), (req, res) => {
    remove((req.params.id as string || '').toLowerCase());
    res.status(200);
    res.end();
});

app.get(mapPath(`/latest`), (req, res) => {
    const index = readIndex();
    const id: string = index.length > 0 ? index[index.length - 1].id : '';
    console.log('Reading', id);
    if (id.length === 0 || id.indexOf('.') >= 0 || id.indexOf('/') >= 0 || id.indexOf('\\') >= 0) {
        res.status(404);
        res.end();
        return;
    }

    fs.readFile(path.join(Config.working_folder, id + '.json'), 'utf-8', (err, data) => {
        if (err) {
            res.status(404);
            res.end();
            return;
        }

        res.writeHead(200, {
            'Content-Type': 'application/json; charset=utf-8',
        });
        res.write(data);
        res.end();
    });
});

app.get(mapPath(`/list`), (req, res) => {
    const index = readIndex();
    res.writeHead(200, {
        'Content-Type': 'application/json; charset=utf-8',
    });
    res.write(JSON.stringify(index, null, 2));
    res.end();
});

app.post(mapPath(`/set`), (req, res) => {
    console.log('SET', req.query.name, req.query.description);
    const index = readIndex();
    validateIndex(index);

    const name = (req.query.name as string || new Date().toUTCString()).substr(0, 50);
    const description = (req.query.description as string || '').substr(0, 100);

    index.push({ timestamp: +new Date(), id: createv4(), name, description });
    const entry = index[index.length - 1];

    const data = JSON.stringify({
        id: entry.id,
        name,
        description,
        data: req.body
    });

    fs.writeFile(path.join(Config.working_folder, entry.id + '.json'), data, { encoding: 'utf8' }, () => res.end());
    writeIndex(index);
});

app.get(`*`, (req, res) => {
    res.writeHead(200, {
        'Content-Type': 'text/plain; charset=utf-8'
    });
    res.write(`
GET /list
GET /get/:id
GET /remove/:id
GET /latest
POST /set?name=...&description=... [JSON data]
`);
    res.end();
})

createIndex();
app.listen(Config.port);

console.log(`Mol* PluginState Server`);
console.log('');
console.log(JSON.stringify(Config, null, 2));