/**
 * Copyright (c) 2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import './index.html';
import { FontAtlas } from '../../mol-geo/geometry/text/font-atlas';
import { printTextureImage } from '../../mol-gl/renderable/util';

function test() {
    console.time('FontAtlas init');
    const fontAtlas = new FontAtlas({ fontQuality: 3 });
    console.timeEnd('FontAtlas init');

    console.time('Basic Latin (subset)');
    for (let i = 0x0020; i <= 0x007E; ++i) fontAtlas.get(String.fromCharCode(i));
    console.timeEnd('Basic Latin (subset)');

    console.time('Latin-1 Supplement (subset)');
    for (let i = 0x00A1; i <= 0x00FF; ++i) fontAtlas.get(String.fromCharCode(i));
    console.timeEnd('Latin-1 Supplement (subset)');

    console.time('Greek and Coptic (subset)');
    for (let i = 0x0391; i <= 0x03C9; ++i) fontAtlas.get(String.fromCharCode(i));
    console.timeEnd('Greek and Coptic (subset)');

    console.time('Cyrillic (subset)');
    for (let i = 0x0400; i <= 0x044F; ++i) fontAtlas.get(String.fromCharCode(i));
    console.timeEnd('Cyrillic (subset)');

    console.time('Angstrom Sign');
    fontAtlas.get(String.fromCharCode(0x212B));
    console.timeEnd('Angstrom Sign');

    printTextureImage(fontAtlas.texture, 0.5);
    console.log(`${Object.keys(fontAtlas.mapped).length} chars prepared`);
}

test();