/**
 * Copyright (c) 2019-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author Gianluca Tomasello <giagitom@gmail.com>
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import { ChunkedArray } from '@molstar/core/data/util';
import { Text } from './text.js';
import { getFontAtlas } from './font-atlas.js';
import { assertUnreachable } from '@molstar/core/util/type-helpers';

// avoiding namespace lookup improved performance in Chrome (Aug 2020)
const caAdd4 = ChunkedArray.add4;

export interface TextBuilder {
  add(str: string, x: number, y: number, z: number, depth: number, scale: number, group: number): void;
  getText(): Text;
}

export namespace TextBuilder {
  export function create(
    props: Partial<PD.Values<Text.Params>> = {},
    initialCount = 2048,
    chunkSize = 1024,
    text?: Text,
  ): TextBuilder {
    // the geometry's data (`TextStride` floats per item: a glyph, background or tether) as 4 elements of 4 floats,
    // reused where given; `initialCount` and `chunkSize` count items
    const data = ChunkedArray.create(
      Float32Array,
      4,
      Math.max(1, chunkSize) * 4,
      text ? text.textData.ref.value.array : Math.max(1, initialCount) * 4,
    );

    const p = { ...PD.getDefaultValues(Text.Params), ...props };
    const { attachment, background, backgroundMargin, tether, tetherLength, tetherBaseWidth } = p;

    const fontAtlas = getFontAtlas(p);
    const margin = (1 / 2.5) * backgroundMargin;
    const outline = fontAtlas.buffer / fontAtlas.lineHeight;

    /** an item of the label at x, y, z: a glyph or a background (`tether` false), or a tether, by its 8 numbers */
    const addItem = (
      x: number,
      y: number,
      z: number,
      depth: number,
      group: number,
      tether: boolean,
      a0: number,
      a1: number,
      a2: number,
      a3: number,
      b0: number,
      b1: number,
      b2: number,
      b3: number,
    ) => {
      caAdd4(data, x, y, z, group);
      caAdd4(data, depth, tether ? 1 : 0, 0, 0);
      caAdd4(data, a0, a1, a2, a3);
      caAdd4(data, b0, b1, b2, b3);
    };

    return {
      add: (str: string, x: number, y: number, z: number, depth: number, scale: number, group: number) => {
        // graphemes placed by measured advances with pair kerning, not by the widths of their atlas cells (labels.md §4)
        const { graphemes, offsets, width: bWidth } = fontAtlas.layout(str);

        const bHeight = 1 / 1.25;

        // attachment
        let yShift: number, xShift: number;
        // vertical
        if (attachment.startsWith('top')) {
          yShift = bHeight;
        } else if (attachment.startsWith('middle')) {
          yShift = bHeight / 2;
        } else {
          yShift = 0; // "bottom"
        }
        // horizontal
        if (attachment.endsWith('right')) {
          xShift = bWidth;
        } else if (attachment.endsWith('center')) {
          xShift = bWidth / 2;
        } else {
          xShift = 0; // "left"
        }

        if (tether) {
          switch (attachment) {
            case 'bottom-left':
              xShift -= tetherLength / 2 + margin + 0.1;
              yShift -= tetherLength / 2 + margin;
              break;
            case 'bottom-center':
              yShift -= tetherLength + margin;
              break;
            case 'bottom-right':
              xShift += tetherLength / 2 + margin + 0.1;
              yShift -= tetherLength / 2 + margin;
              break;
            case 'middle-left':
              xShift -= tetherLength + margin + 0.1;
              break;
            case 'middle-center':
              break;
            case 'middle-right':
              xShift += tetherLength + margin + 0.1;
              break;
            case 'top-left':
              xShift -= tetherLength / 2 + margin + 0.1;
              yShift += tetherLength / 2 + margin;
              break;
            case 'top-center':
              yShift += tetherLength + margin;
              break;
            case 'top-right':
              xShift += tetherLength / 2 + margin + 0.1;
              yShift += tetherLength / 2 + margin;
              break;
          }
        }

        const xLeft = (-xShift - margin - 0.1) * scale;
        const xRight = (bWidth - xShift + margin + 0.1) * scale;
        const yMid = 0.5 - yShift - outline * 0.5; // glyph vertical midpoint accounting for outline offset
        const yTop = (yMid + bHeight / 2 + margin) * scale;
        const yBottom = (yMid - bHeight / 2 - margin) * scale;

        // background: texture coordinates 10 (no glyph)
        if (background) {
          addItem(x, y, z, depth, group, false, xLeft, yTop, xRight, yBottom, 10, 10, 10, 10);
        }

        if (tether) {
          let xTip: number, yTip: number;
          let xBaseA: number, yBaseA: number;
          let xBaseB: number, yBaseB: number;
          let xBaseCenter: number, yBaseCenter: number;

          const scaledTetherLength = tetherLength * scale;
          const scaledTetherBaseWidth = tetherBaseWidth * scale;

          switch (attachment) {
            case 'bottom-left':
              xTip = xLeft - scaledTetherLength / 2;
              xBaseA = xLeft + scaledTetherBaseWidth / 2;
              xBaseB = xLeft;
              xBaseCenter = xLeft;
              yTip = yBottom - scaledTetherLength / 2;
              yBaseA = yBottom;
              yBaseB = yBottom + scaledTetherBaseWidth / 2;
              yBaseCenter = yBottom;
              break;
            case 'bottom-center':
              xTip = 0;
              xBaseA = scaledTetherBaseWidth / 2;
              xBaseB = -scaledTetherBaseWidth / 2;
              xBaseCenter = 0;
              yTip = yBottom - scaledTetherLength;
              yBaseA = yBottom;
              yBaseB = yBottom;
              yBaseCenter = yBottom;
              break;
            case 'bottom-right':
              xTip = xRight + scaledTetherLength / 2;
              xBaseA = xRight;
              xBaseB = xRight - scaledTetherBaseWidth / 2;
              xBaseCenter = xRight;
              yTip = yBottom - scaledTetherLength / 2;
              yBaseA = yBottom + scaledTetherBaseWidth / 2;
              yBaseB = yBottom;
              yBaseCenter = yBottom;
              break;
            case 'middle-left':
              xTip = xLeft - scaledTetherLength;
              xBaseA = xLeft;
              xBaseB = xLeft;
              xBaseCenter = xLeft;
              yTip = 0;
              yBaseA = -scaledTetherBaseWidth / 2;
              yBaseB = scaledTetherBaseWidth / 2;
              yBaseCenter = 0;
              break;
            case 'middle-center':
              xTip = 0;
              xBaseA = 0;
              xBaseB = 0;
              xBaseCenter = 0;
              yTip = 0;
              yBaseA = 0;
              yBaseB = 0;
              yBaseCenter = 0;
              break;
            case 'middle-right':
              xTip = xRight + scaledTetherLength;
              xBaseA = xRight;
              xBaseB = xRight;
              xBaseCenter = xRight;
              yTip = 0;
              yBaseA = scaledTetherBaseWidth / 2;
              yBaseB = -scaledTetherBaseWidth / 2;
              yBaseCenter = 0;
              break;
            case 'top-left':
              xTip = xLeft - scaledTetherLength / 2;
              xBaseA = xLeft + scaledTetherBaseWidth / 2;
              xBaseB = xLeft;
              xBaseCenter = xLeft;
              yTip = yTop + scaledTetherLength / 2;
              yBaseA = yTop;
              yBaseB = yTop - scaledTetherBaseWidth / 2;
              yBaseCenter = yTop;
              break;
            case 'top-center':
              xTip = 0;
              xBaseA = scaledTetherBaseWidth / 2;
              xBaseB = -scaledTetherBaseWidth / 2;
              xBaseCenter = 0;
              yTip = yTop + scaledTetherLength;
              yBaseA = yTop;
              yBaseB = yTop;
              yBaseCenter = yTop;
              break;
            case 'top-right':
              xTip = xRight + scaledTetherLength / 2;
              xBaseA = xRight;
              xBaseB = xRight - scaledTetherBaseWidth / 2;
              xBaseCenter = xRight;
              yTip = yTop + scaledTetherLength / 2;
              yBaseA = yTop - scaledTetherBaseWidth / 2;
              yBaseB = yTop;
              yBaseCenter = yTop;
              break;
            default:
              assertUnreachable(attachment);
          }
          // the corners in the order the shared indices draw the triangles (tip, base A, base center) and (tip,
          // base center, base B) in: base B, tip, base center, base A
          addItem(x, y, z, depth, group, true, xBaseB, yBaseB, xTip, yTip, xBaseCenter, yBaseCenter, xBaseA, yBaseA);
        }

        yShift += outline;

        for (let iChar = 0, nChar = graphemes.length; iChar < nChar; ++iChar) {
          const c = fontAtlas.get(graphemes[iChar]);
          const left = (offsets[iChar] - c.left - xShift) * scale;
          const right = left + c.nw * scale;
          const top = (c.nh - yShift) * scale;
          const bottom = -yShift * scale;

          const texWidth = fontAtlas.texture.width;
          const texHeight = fontAtlas.texture.height;

          // corners and texture coordinates as left, top, right, bottom
          addItem(
            x,
            y,
            z,
            depth,
            group,
            false,
            left,
            top,
            right,
            bottom,
            c.x / texWidth,
            c.y / texHeight,
            (c.x + c.w) / texWidth,
            (c.y + c.h) / texHeight,
          );
        }
      },
      getText: () => {
        const array = ChunkedArray.compactTrimmed(data) as Float32Array;
        return Text.create(fontAtlas.texture, array, data.elementCount / 4, text);
      },
    };
  }
}
