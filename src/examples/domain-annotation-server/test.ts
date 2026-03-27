/**
 * Copyright (c) 2017-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Paul Pillot <paul.pillot@tandemai.com>
 */

import { createMapping } from './mapping';

(async function () {
    const data = await fetch('https://www.ebi.ac.uk/pdbe/api/mappings/1tqn?pretty=true');
    const json = await data.json();
    console.log(createMapping(json));
}());