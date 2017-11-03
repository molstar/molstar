/**
 * Copyright (c) 2017 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Adapted from CIFTools.js (https://github.com/dsehnal/CIFTools.js)
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */


interface StringBuilder {
    current: string[],
    offset: number,
    capacity: number,
    chunks: string[]
}

namespace StringBuilder {
    export function create(chunkCapacity = 512): StringBuilder {
        return {
            current: [],
            offset: 0,
            capacity: chunkCapacity,
            chunks: []
        };
    }

    export function getString(builder: StringBuilder) {
        if (!builder.chunks.length) {
            if (builder.current.length === builder.offset) return builder.current.join('');
            return builder.current.splice(0, builder.offset).join('');
        }

        if (builder.offset > 0) {
            builder.chunks[builder.chunks.length] = builder.current.length === builder.offset
                ? builder.current.join('')
                : builder.current.slice(0, builder.offset).join('');
        }

        return builder.chunks.join('');
    }

    export function getChunks(builder: StringBuilder): string[] {
        if (builder.offset > 0) {
            if (builder.current.length === builder.offset) builder.chunks[builder.chunks.length] = builder.current.join('');
            else builder.chunks[builder.chunks.length] = builder.current.slice(0, builder.offset).join('');
            builder.offset = 0;
        }
        return builder.chunks;
    }

    const enum PaddingSpaces { Count = 512 }
    const __paddingSpaces: string[] = [];
    (function () {
        let s = '';
        for (let i = 0; i < PaddingSpaces.Count; i++) {
            __paddingSpaces[i] = s;
            s = s + ' ';
        }
    })();

    export function newline(builder: StringBuilder) {
        writeSafe(builder, '\n');
    }

    export function whitespace(builder: StringBuilder, len: number) {
        if (len > 0) write(builder, __paddingSpaces[len]);
    }

    export function write(builder: StringBuilder, val: string) {
        if (!val) return;

        if (builder.offset === builder.capacity) {
            builder.chunks[builder.chunks.length] = builder.current.join('');
            builder.offset = 0;
        }

        builder.current[builder.offset++] = val;
    }

    /** Write without check. */
    export function writeSafe(builder: StringBuilder, val: string) {
        if (builder.offset === builder.capacity) {
            builder.chunks[builder.chunks.length] = builder.current.join('');
            builder.offset = 0;
        }

        builder.current[builder.offset++] = val;
    }

    export function writePadLeft(builder: StringBuilder, val: string, totalWidth: number) {
        if (!val) { whitespace(builder, totalWidth); return; }

        let padding = totalWidth - val.length;
        whitespace(builder, padding);
        writeSafe(builder, val);
    }

    export function writePadRight(builder: StringBuilder, val: string, totalWidth: number) {
        if (!val) { whitespace(builder, totalWidth); return; }

        let padding = totalWidth - val.length;
        writeSafe(builder, val);
        whitespace(builder, padding);
    }


    export function writeInteger(builder: StringBuilder, val: number) {
        writeSafe(builder, '' + val);
    }

    export function writeIntegerPadLeft(builder: StringBuilder, val: number, totalWidth: number) {
        let s = '' + val;
        let padding = totalWidth - s.length;
        whitespace(builder, padding);
        writeSafe(builder, s);
    }

    export function writeIntegerPadRight(builder: StringBuilder, val: number, totalWidth: number) {
        let s = '' + val;
        let padding = totalWidth - s.length;
        writeSafe(builder, s);
        whitespace(builder, padding);
    }

    /**
     * @example writeFloat(123.2123, 100) -- 2 decim
     */
    export function writeFloat(builder: StringBuilder, val: number, precisionMultiplier: number) {
        writeSafe(builder, '' + Math.round(precisionMultiplier * val) / precisionMultiplier)
    }

    export function writeFloatPadLeft(builder: StringBuilder, val: number, precisionMultiplier: number, totalWidth: number) {
        let s = '' + Math.round(precisionMultiplier * val) / precisionMultiplier;
        let padding = totalWidth - s.length;
        whitespace(builder, padding);
        writeSafe(builder, s);
    }

    export function writeFloatPadRight(builder: StringBuilder, val: number, precisionMultiplier: number, totalWidth: number) {
        let s = '' + Math.round(precisionMultiplier * val) / precisionMultiplier;
        let padding = totalWidth - s.length;
        writeSafe(builder, s);
        whitespace(builder, padding);
    }
}

export default StringBuilder