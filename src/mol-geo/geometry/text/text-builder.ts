/**
 * Copyright (c) 2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { ParamDefinition as PD } from 'mol-util/param-definition';
import { ValueCell } from 'mol-util/value-cell'
import { ChunkedArray } from 'mol-data/util';
import { Text } from './text';
import { getFontAtlas } from './font-atlas';

const quadIndices = new Uint16Array([
    0, 1, 2,
    1, 3, 2
])

export interface TextBuilder {
    add(str: string, x: number, y: number, z: number, depth: number, group: number): void
    getText(): Text
}

export namespace TextBuilder {
    export function create(props: Partial<PD.Values<Text.Params>> = {}, initialCount = 2048, chunkSize = 1024, text?: Text): TextBuilder {
        initialCount *= 2
        chunkSize *= 2
        const centers = ChunkedArray.create(Float32Array, 3, chunkSize, text ? text.centerBuffer.ref.value : initialCount);
        const mappings = ChunkedArray.create(Float32Array, 2, chunkSize, text ? text.mappingBuffer.ref.value : initialCount);
        const depths = ChunkedArray.create(Float32Array, 1, chunkSize, text ? text.depthBuffer.ref.value : initialCount);
        const indices = ChunkedArray.create(Uint32Array, 3, chunkSize, text ? text.indexBuffer.ref.value : initialCount);
        const groups = ChunkedArray.create(Float32Array, 1, chunkSize, text ? text.groupBuffer.ref.value : initialCount);
        const tcoords = ChunkedArray.create(Float32Array, 2, chunkSize, text ? text.tcoordBuffer.ref.value : initialCount);

        const p = { ...PD.getDefaultValues(Text.Params), ...props }
        const { attachment, background, backgroundMargin } = p

        const fontAtlas = getFontAtlas(p)
        const margin = (1 / 2.5) * backgroundMargin
        const outline = fontAtlas.buffer / fontAtlas.lineHeight

        return {
            add: (str: string, x: number, y: number, z: number, depth: number, group: number) => {
                let bWidth = 0
                const nChar = str.length

                // calculate width
                for (let iChar = 0; iChar < nChar; ++iChar) {
                    const c = fontAtlas.get(str[iChar])
                    bWidth += c.nw - 2 * outline
                }

                const bHeight = 1 / 1.25

                // attachment
                let yShift: number, xShift: number
                if (attachment.startsWith('top')) {
                    yShift = bHeight
                } else if (attachment.startsWith('middle')) {
                    yShift = bHeight / 2
                } else {
                    yShift = 0  // "bottom"
                }
                if (attachment.endsWith('right')) {
                    xShift = bWidth
                } else if (attachment.endsWith('center')) {
                    xShift = bWidth / 2
                } else {
                    xShift = 0  // "left"
                }

                // background
                if (background) {
                    ChunkedArray.add2(mappings, -xShift - margin - 0.1, bHeight - yShift + margin) // top left
                    ChunkedArray.add2(mappings, -xShift - margin - 0.1, -yShift - margin) // bottom left
                    ChunkedArray.add2(mappings, bWidth - xShift + margin + 0.1, bHeight - yShift + margin) // top right
                    ChunkedArray.add2(mappings, bWidth - xShift + margin + 0.1, -yShift - margin) // bottom right

                    const offset = centers.elementCount
                    for (let i = 0; i < 4; ++i) {
                        ChunkedArray.add2(tcoords, 10, 10)
                        ChunkedArray.add3(centers, x, y, z);
                        ChunkedArray.add(depths, depth);
                        ChunkedArray.add(groups, group);
                    }
                    ChunkedArray.add3(indices, offset + quadIndices[0], offset + quadIndices[1], offset + quadIndices[2])
                    ChunkedArray.add3(indices, offset + quadIndices[3], offset + quadIndices[4], offset + quadIndices[5])
                }

                xShift += outline
                yShift += outline
                let xadvance = 0

                for (let iChar = 0; iChar < nChar; ++iChar) {
                    const c = fontAtlas.get(str[iChar])

                    ChunkedArray.add2(mappings, xadvance - xShift, c.nh - yShift) // top left
                    ChunkedArray.add2(mappings, xadvance - xShift, -yShift) // bottom left
                    ChunkedArray.add2(mappings, xadvance + c.nw - xShift, c.nh - yShift) // top right
                    ChunkedArray.add2(mappings, xadvance + c.nw - xShift, -yShift) // bottom right

                    const texWidth = fontAtlas.texture.width
                    const texHeight = fontAtlas.texture.height

                    ChunkedArray.add2(tcoords, c.x / texWidth, c.y / texHeight) // top left
                    ChunkedArray.add2(tcoords, c.x / texWidth, (c.y + c.h) / texHeight) // bottom left
                    ChunkedArray.add2(tcoords, (c.x + c.w) / texWidth, c.y / texHeight) // top right
                    ChunkedArray.add2(tcoords, (c.x + c.w) / texWidth, (c.y + c.h) / texHeight) // bottom right

                    xadvance += c.nw - 2 * outline

                    const offset = centers.elementCount
                    for (let i = 0; i < 4; ++i) {
                        ChunkedArray.add3(centers, x, y, z);
                        ChunkedArray.add(depths, depth);
                        ChunkedArray.add(groups, group);
                    }
                    ChunkedArray.add3(indices, offset + quadIndices[0], offset + quadIndices[1], offset + quadIndices[2])
                    ChunkedArray.add3(indices, offset + quadIndices[3], offset + quadIndices[4], offset + quadIndices[5])
                }
            },
            getText: () => {
                const ft = fontAtlas.texture
                const cb = ChunkedArray.compact(centers, true) as Float32Array
                const mb = ChunkedArray.compact(mappings, true) as Float32Array
                const db = ChunkedArray.compact(depths, true) as Float32Array
                const ib = ChunkedArray.compact(indices, true) as Uint32Array
                const gb = ChunkedArray.compact(groups, true) as Float32Array
                const tb = ChunkedArray.compact(tcoords, true) as Float32Array
                return {
                    kind: 'text',
                    charCount: indices.elementCount / 2,
                    fontTexture: text ? ValueCell.update(text.fontTexture, ft) : ValueCell.create(ft),
                    centerBuffer: text ? ValueCell.update(text.centerBuffer, cb) : ValueCell.create(cb),
                    mappingBuffer: text ? ValueCell.update(text.mappingBuffer, mb) : ValueCell.create(mb),
                    depthBuffer: text ? ValueCell.update(text.depthBuffer, db) : ValueCell.create(db),
                    indexBuffer: text ? ValueCell.update(text.indexBuffer, ib) : ValueCell.create(ib),
                    groupBuffer: text ? ValueCell.update(text.groupBuffer, gb) : ValueCell.create(gb),
                    tcoordBuffer: text ? ValueCell.update(text.tcoordBuffer, tb) : ValueCell.create(tb),
                }
            }
        }
    }
}