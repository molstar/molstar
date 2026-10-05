/**
 * Copyright (c) 2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import './index.html';
import { FontAtlas } from '@molstar/graphics/geo/geometry/text/font-atlas';
import { printTextureImage } from '@molstar/graphics/gl/renderable/util';

function test() {
  console.time('FontAtlas init');
  const fontAtlas = new FontAtlas({ fontQuality: 3 });
  console.timeEnd('FontAtlas init');

  console.time('Basic Latin (subset)');
  for (let i = 0x0020; i <= 0x007e; ++i) fontAtlas.get(String.fromCharCode(i));
  console.timeEnd('Basic Latin (subset)');

  console.time('Latin-1 Supplement (subset)');
  for (let i = 0x00a1; i <= 0x00ff; ++i) fontAtlas.get(String.fromCharCode(i));
  console.timeEnd('Latin-1 Supplement (subset)');

  console.time('Greek and Coptic (subset)');
  for (let i = 0x0391; i <= 0x03c9; ++i) fontAtlas.get(String.fromCharCode(i));
  console.timeEnd('Greek and Coptic (subset)');

  console.time('Cyrillic (subset)');
  for (let i = 0x0400; i <= 0x044f; ++i) fontAtlas.get(String.fromCharCode(i));
  console.timeEnd('Cyrillic (subset)');

  console.time('Angstrom Sign');
  fontAtlas.get(String.fromCharCode(0x212b));
  console.timeEnd('Angstrom Sign');

  printTextureImage(fontAtlas.texture, { scale: 0.5 });
  console.log(`${Object.keys(fontAtlas.mapped).length} chars prepared`);
}

test();
