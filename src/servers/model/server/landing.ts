/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import Version from '../version'

const examples = [{
    name: 'Atoms',
    params: {
        id: '1cbs',
        name: 'atoms',
        params: { atom_site: { label_comp_id: 'ALA' } }
    }
}, {
    name: 'Residue Interaction',
    params: {
        id: '1cbs',
        name: 'residueInteraction',
        params: {
            radius: 5,
            atom_site: { 'label_comp_id': 'REA' }
        }
    }
}, {
    name: 'Full',
    params: {
        id: '1tqn',
        name: 'full'
    }
}, {
    name: 'Full (binary)',
    params: {
        id: '1tqn',
        name: 'full',
        binary: true
    }
}, {
    name: 'Full (specific models)',
    params: {
        id: '1grm',
        name: 'full',
        modelNums: [ 2, 3 ]
    }
}];

function create() {
    return `<!DOCTYPE html>
<html lang="en">
    <head>
        <meta charset="utf-8" />
        <meta name="viewport" content="width=device-width, user-scalable=no, minimum-scale=1.0, maximum-scale=1.0">
        <title>Mol* ModelServer ${Version}</title>
        <link rel="stylesheet" href="https://cdnjs.cloudflare.com/ajax/libs/skeleton/2.0.4/skeleton.min.css" />
    </head>
    <body>
        <h1>Mol* Model Server ${Version}</h1>
        <select id='example'>
            <option value='-1'>Select example...</option>
            ${examples.map((e, i) => `<option value=${i}>${e.name}</option>`)}
        </select>
        <br/>
        <textarea style="height: 280px; width: 600px; font-family: monospace" id="query-text"></textarea><br>
        <button class="button button-primary" style="width: 600px" id="query">Query</button>
        <div id='error' style='color: red; font-weight: blue'></div>
        <div>Static input files available as CIF and BinaryCIF at <a href='ModelServer/static/cif/1cbs' target='_blank'>static/cif/id</a> and <a href='ModelServer/static/bcif/1cbs' target='_blank'>static/bcif/id</a> respectively.</div>
        <script>
            var Examples = ${JSON.stringify(examples)};
            var err = document.getElementById('error');
            var exampleEl = document.getElementById('example'), queryTextEl = document.getElementById('query-text');
            exampleEl.onchange = function () {
                var i = +exampleEl.value;
                if (i < 0) return;
                queryTextEl.value = JSON.stringify(Examples[i].params, null, 2);
            };
            document.getElementById('query').onclick = function () {
                err.innerText = '';
                try {
                    var q = JSON.parse(queryTextEl.value);
                    var path = 'ModelServer/api/v1?' + encodeURIComponent(JSON.stringify(q));
                    console.log(path);
                    window.open(path, '_blank');
                } catch (e) {
                    err.innerText = '' + e;
                }
            };
        </script>
    </body>
</html>`;
}

export const LandingPage = create();