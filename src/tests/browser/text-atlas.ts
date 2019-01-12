/**
 * Copyright (c) 2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import './index.html'
import { TextAtlas } from 'mol-geo/geometry/text/text-atlas';
import { printImageData } from 'mol-gl/renderable/util';

function test() {
    console.time('TextAtlas')
    const textAtlas = new TextAtlas()
    console.timeEnd('TextAtlas')
    const ctx = textAtlas.context
    const imageData = ctx.getImageData(0, 0, ctx.canvas.width, ctx.canvas.height)
    printImageData(imageData)
}

test();