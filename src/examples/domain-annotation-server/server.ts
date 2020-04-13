/**
 * Copyright (c) 2017 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import express from 'express';
import fetch from 'node-fetch';
import createMapping from './mapping';

async function getMappings(id: string) {
    const data = await fetch(`https://www.ebi.ac.uk/pdbe/api/mappings/${id}`);
    const json = await data.json();
    return createMapping(json);
};


let PORT = process.env.port || 1338;

const app = express();

const PREFIX = '/';

app.get(`${PREFIX}/:id`, async (req, res) => {
    try {
        console.log('Requesting ' + req.params.id);
        const mapping = await getMappings(req.params.id);
        res.writeHead(200, {
            'Content-Type': 'text/plain; charset=utf-8',
            'Access-Control-Allow-Origin': '*',
            'Access-Control-Allow-Headers': 'X-Requested-With'
        });
        res.end(mapping);
    } catch {
        console.log('Failed ' + req.params.id);
        res.writeHead(404, { 'Access-Control-Allow-Origin': '*', 'Access-Control-Allow-Headers': 'X-Requested-With' });
        res.end();
    }
});

app.get(`${PREFIX}`, (req, res) => {
    res.writeHead(200, { 'Content-Type': 'text/plain; charset=utf-8' });
    res.end('Usage: /pdb_id, e.g. /1tqn');
});

app.listen(PORT);

console.log('Running on port ' + PORT);