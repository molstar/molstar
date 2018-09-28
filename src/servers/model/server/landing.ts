/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import Version from '../version'

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
        <textarea style="height: 280px; width: 600px; font-family: monospace" id="query-text">{
    "id": "1cbs",
    "name": "residueInteraction",
    "params": {
        "radius": 5,
        "atom_site": { "label_comp_id": "REA" }
    }
}</textarea><br>
        <button class="button button-primary" style="width: 600px" id="query">Query</button>
        <div id='error' style='color: red; font-weight: blue'></div>
        <div>Static input files available as CIF and BinaryCIF at <a href='ModelServer/static/cif/1cbs' target='_blank'>static/cif/id</a> and <a href='ModelServer/static/bcif/1cbs' target='_blank'>static/bcif/id</a> respectively.</div>
        <script>
            const err = document.getElementById('error');
            document.getElementById('query').onclick = function () {
                err.innerText = '';
                try {
                    var q = JSON.parse(document.getElementById('query-text').value);
                    var path = 'ModelServer/api/v1?' + encodeURIComponent(JSON.stringify(q));
                    console.log(path);
                    window.open(path, '_blank');
                } catch (e) {
                    err.innerText = '' + e;
                }
            }
        </script>
    </body>
</html>`;
}

export const LandingPage = create();