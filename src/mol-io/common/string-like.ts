/**
 * Copyright (c) 2025 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Adam Midlik <midlik@gmail.com>
 */

import { utf8Read } from './utf8';


/**
 * Essential subset of `string` functionality.
 * Can be builtin `string` or `String` type or a class instance implementing necessary methods (`StringLikeInterface` interface).
 * Add more string methods if needed.
 */
export type StringLike = string | String | StringLikeInterface; // using this type union because with some TypeScript versions `string` cannot be assigned to `StringLikeInterface`

/** Classes that want to implement `StringLike` should use `implements StringLikeInterface` (it is not possible to use `implements` with union type directly). */
export interface StringLikeInterface {
    /** Returns the length of a String object. */
    readonly length: number;

    // readonly [index: number]: string; // Avoided implementing this as it would require using Proxy, which drastically affects performance.

    /**
     * Returns the character at the specified index, or `undefined` if the index is out of range. Supports relative indexing from the end of the string when passed a negative index.
     * @param pos The zero-based index of the desired character.
     */
    at(index: number): string | undefined;

    /**
     * Returns the character at the specified index, or an empty string if the index is out of range.
     * @param pos The zero-based index of the desired character.
     */
    charAt(pos: number): string;

    /**
     * Returns the Unicode value of the character at the specified location.
     * @param index The zero-based index of the desired character. If the specified index is out of range, NaN is returned.
     */
    charCodeAt(index: number): number;

    /**
     * Returns the substring at the specified location within a String object.
     * @param start The zero-based index number indicating the beginning of the substring.
     * @param end Zero-based index number indicating the end of the substring. The substring includes the characters up to, but not including, the character indicated by end.
     * If end is omitted, the characters from start through the end of the original string are returned.
     */
    substring(start: number, end?: number): string;

    /**
     * Returns the position of the first occurrence of a substring, or -1 if not found.
     * @param searchString The substring to search for in the string
     * @param position The index at which to begin searching the String object. If omitted, search starts at the beginning of the string.
     */
    indexOf(searchString: string, position?: number): number;

    /**
     * Returns true if searchString appears as a substring of the result of converting this
     * object to a String, at one or more positions that are
     * greater than or equal to position; otherwise, returns false.
     * @param searchString search string
     * @param position If position is undefined, 0 is assumed, so as to search all of the String.
     */
    includes(searchString: string, position?: number): boolean;

    /**
     * Returns true if the sequence of elements of searchString converted to a String is the
     * same as the corresponding elements of this object (converted to a String) starting at
     * position. Otherwise returns false.
     */
    startsWith(searchString: string, position?: number): boolean;

    /** Returns a string representation of a string. */
    toString(): string;
}

export const StringLike = {
    /** Return true if `obj` is instance of `StringLike` */
    is(obj: unknown): obj is StringLike {
        return typeof (obj as StringLike).charCodeAt === 'function'; // a bit hacky
    },

    /** Try to convert `StringLike` to a primitive `string`. Might fail if the content is longer that max allowed string length. */
    toString(str: StringLike): string {
        try {
            return str.toString();
        } catch (err) {
            throw new Error(`Failed to convert StringLike object into string. This might be because the length ${str.length} exceeds maximum allowed string length ${MAX_STRING_LENGTH}. (${err})`);
        }
    },
};


/** Maximum allowed string length (might be bigger for some engines, but in Chrome 136 and Node 22 it is this). */
export const MAX_STRING_LENGTH = 536_870_888;

/** Binary logarithm of default string chunk size for `ChunkedBigString`. (string chunk size is chosen to be a power of 2, so we can use faster bit shift operator instead of integer division) */
const DEFAULT_LOG_STRING_CHUNK_SIZE = 28; // 2**28 is the largest power of 2 which is <= MAX_STRING_LENGTH


/** Implementation of `CustomString`, based on an array of fixed-length strings (chunks). */
export class ChunkedBigString implements StringLikeInterface {
    private _chunks: string[] = [];

    /** Length of string chunks (default 2**28). */
    private readonly STRING_CHUNK_SIZE: number;
    /** Bit shift to get chunk index from char index (default 28) */
    private readonly STRING_CHUNK_SHIFT: number;
    /** Bit mask to get index within chunk index from char index (default 2**28 - 1) */
    private readonly STRING_CHUNK_MASK: number;

    private _length: number = 0;

    get length(): number {
        return this._length;
    }

    constructor(logStringChunkSize: number = DEFAULT_LOG_STRING_CHUNK_SIZE) {
        this.STRING_CHUNK_SIZE = 2 ** logStringChunkSize;
        this.STRING_CHUNK_SHIFT = logStringChunkSize;
        this.STRING_CHUNK_MASK = 2 ** logStringChunkSize - 1;
    }

    static fromString(content: string, logStringChunkSize: number = DEFAULT_LOG_STRING_CHUNK_SIZE): ChunkedBigString {
        const out = new ChunkedBigString(logStringChunkSize);
        out._append(content);
        return out;
    }

    static fromStrings(content: string[], logStringChunkSize: number = DEFAULT_LOG_STRING_CHUNK_SIZE): ChunkedBigString {
        const out = new ChunkedBigString(logStringChunkSize);
        for (const inputChunk of content) {
            out._append(inputChunk);
        }
        return out;
    }

    /** Create instance from UTF8 data. (Do not call directly, prefer `utf8ReadLong` in utf8.ts.) */
    static fromUtf8Data(data: Uint8Array, start: number = 0, end: number = data.length, logStringChunkSize: number = DEFAULT_LOG_STRING_CHUNK_SIZE): ChunkedBigString {
        const bufferChunkSize = 2 ** logStringChunkSize; // n bytes will always decode to <=n characters
        const stringChunks: string[] = [];
        let readStart = start;
        while (readStart < end) {
            let readEnd = Math.min(readStart + bufferChunkSize, end);
            if (readEnd < end) {
                // This is buffer chunk boundary, adjust to avoid cutting multi-byte characters
                while ((data[readEnd] & 0xC0) === 0x80) { // Byte after the cut is a continuation byte (10xxxxxx)
                    readEnd--;
                    if (readEnd === readStart) throw new Error('Input is rubbish, no UTF-8 character start found in a chunk');
                }
            } // Else this is the end of the read region, let default error handling do its job
            const stringChunk = utf8Read(data, readStart, readEnd - readStart);
            stringChunks.push(stringChunk);
            readStart = readEnd;
        }
        return ChunkedBigString.fromStrings(stringChunks, logStringChunkSize);
    }

    private _append(inputChunk: string): void {
        const chunkSize = this.STRING_CHUNK_SIZE;
        const tail = (this._chunks.length === 0 || this._chunks[this._chunks.length - 1].length === chunkSize) ? '' : this._chunks.pop()!;
        let inputPtr = chunkSize - tail.length;
        this._chunks.push(tail + inputChunk.substring(0, inputPtr)); // Assuming .substring() deals with inputPtr > inputChunk.length
        while (inputPtr < inputChunk.length) {
            this._chunks.push(inputChunk.substring(inputPtr, inputPtr + chunkSize)); // Assuming .substring() deals with inputPtr + chunkSize > inputChunk.length
            inputPtr += chunkSize;
        }
        this._length += inputChunk.length;
    }

    private _getChunkIndex(index: number) {
        return index >>> this.STRING_CHUNK_SHIFT; // equivalent to `Math.floor(index / STRING_CHUNK_SIZE)`
    }
    private _getIndexInChunk(index: number) {
        return index & this.STRING_CHUNK_MASK; // equivalent to `index % STRING_CHUNK_SIZE`
    }
    private _isOutOfRange(index: number) {
        return index < 0 || index >= this.length;
    }

    at(index: number): string | undefined {
        if (-this.length <= index && index < 0) {
            return this.at(index + this.length);
        }
        return this.charAt(index) || undefined;
    }

    charAt(index: number): string {
        if (this._isOutOfRange(index)) return '';
        const iChunk = this._getChunkIndex(index);
        const indexInChunk = this._getIndexInChunk(index);
        return this._chunks[iChunk][indexInChunk];
    }

    charCodeAt(index: number): number {
        if (this._isOutOfRange(index)) return NaN;
        const iChunk = this._getChunkIndex(index);
        const indexInChunk = this._getIndexInChunk(index);
        return this._chunks[iChunk].charCodeAt(indexInChunk);
    }

    substring(start?: number, end?: number): string { // optional `start` not part of contract but works in Chrome
        const start_ = Math.min(Math.max(start ?? 0, 0), this.length);
        const end_ = Math.min(Math.max(end ?? this.length, 0), this.length);
        if (start_ > end_) {
            return this.substring(end_, start_);
        }
        if (start_ === end_) {
            return '';
        }

        if (end_ - start_ > MAX_STRING_LENGTH) {
            throw new Error(`Trying to create get a substring longer (${end_ - start_}) than maximum allowed string length (${MAX_STRING_LENGTH}).`);
        }

        const iFirstChunk = this._getChunkIndex(start_);
        const indexInChunkFrom = this._getIndexInChunk(start_);
        const iLastChunk = this._getChunkIndex(end_);
        const indexInChunkTo = this._getIndexInChunk(end_);

        if (iFirstChunk === iLastChunk) {
            return this._chunks[iFirstChunk].substring(indexInChunkFrom, indexInChunkTo);
        } else {
            const out = this._getTmpArray();
            out.push(this._chunks[iFirstChunk].substring(indexInChunkFrom, this.STRING_CHUNK_SIZE));
            for (let iChunk = iFirstChunk + 1; iChunk < iLastChunk; iChunk++) {
                out.push(this._chunks[iChunk]);
            }
            out.push(this._chunks[iLastChunk].substring(0, indexInChunkTo));
            return out.join('');
        }
    }
    private readonly _tmpArray: string[] = [];
    private _getTmpArray() {
        while (this._tmpArray.length) this._tmpArray.pop(); // this seems to be faster than `this._tmpArray.length = 0` for short arrays
        return this._tmpArray;
    }

    indexOf(searchString: string, position: number = 0): number {
        if (searchString.length > this.STRING_CHUNK_SIZE) {
            throw new Error('NotImplementedError: indexOf is only implemented for searchString shorter than STRING_CHUNK_SIZE');
            // In real use-cases STRING_CHUNK_SIZE is big and it doesn't make sense to search for such long substrings.
        }

        if (position < 0) position = 0;
        const iFirstChunk = this._getChunkIndex(position);

        for (let iChunk = iFirstChunk; iChunk < this._chunks.length; iChunk++) {
            const chunk = this._chunks[iChunk];
            const positionInChunk = iChunk === iFirstChunk ? this._getIndexInChunk(position) : 0;

            // Try to find the whole substring in this chunk
            const found = chunk.indexOf(searchString, positionInChunk);
            if (found >= 0) return iChunk * this.STRING_CHUNK_SIZE + found;

            // Try to find the substring overflowing to the next chunk (assumes searchString.length <= STRING_CHUNK_SIZE)
            if (iChunk !== this._chunks.length - 1) {
                const start = Math.max(this.STRING_CHUNK_SIZE - searchString.length + 1, positionInChunk);
                const aroundBoundary = chunk.substring(start, undefined) + this._chunks[iChunk + 1].substring(0, searchString.length - 1);
                const found = aroundBoundary.indexOf(searchString);
                if (found >= 0) return iChunk * this.STRING_CHUNK_SIZE + start + found;
            }
        }
        return -1;
    }

    includes(searchString: string, position: number = 0): boolean {
        return this.indexOf(searchString, position) >= 0;
    }

    startsWith(searchString: string, position: number = 0): boolean {
        if (searchString.length > this.STRING_CHUNK_SIZE) {
            throw new Error('NotImplementedError: startsWith is only implemented for searchString shorter than STRING_CHUNK_SIZE');
            // In real use-cases STRING_CHUNK_SIZE is big and it doesn't make sense to search for such long substrings.
        }
        return this.substring(position, position + searchString.length) === searchString;
    }

    toString(): string {
        try {
            return this._chunks.join('');
        } catch (err) {
            throw new Error(`Failed to convert StringLike object into string. This might be because the length ${this.length} exceeds maximum allowed string length ${MAX_STRING_LENGTH}. (${err})`);
        }
    }
}
