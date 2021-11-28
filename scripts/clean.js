/**
 * Copyright (c) 2021 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Michal Malý <malym@ibt.cas.cz>
 */

const fs = require('fs');
const path = require('path');

const toClean = [
    path.resolve(__dirname, '../build'),
    path.resolve(__dirname, '../lib'),
    path.resolve(__dirname, '../tsconfig.tsbuildinfo'),
];

toClean.forEach(ph => {
    if (fs.existsSync(ph))
        fs.rm(ph, { recursive: true }, err => { if (err) console.warn(`Failed to delete ${ph}: ${err}`); });
});
