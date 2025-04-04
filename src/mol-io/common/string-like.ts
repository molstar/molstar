/**
 * Copyright (c) 2025 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Adam Midlik <midlik@gmail.com>
 */


interface CustomString {
    // TODO think about this extending full `String` interface, and impl classes just throwing `NotImplementedError` on more advanced methods -> less likely to break typing of dependant code (will only break runtime when running on big strings)
    __string_like__: true,

    readonly [index: number]: string;

    // [Symbol.iterator]: () => StringIterator<string>;

    // from @types/node/compatibility/indexable.d.ts

    at(index: number): string | undefined;

    // // from typescript/lib/lib.es5.d.ts

    /** Returns a string representation of a string. */
    toString(): string; // TODO make supersure this is implemented properly

    // /**
    //  * Returns the character at the specified index.
    //  * @param pos The zero-based index of the desired character.
    //  */
    // charAt(pos: number): string;

    /**
     * Returns the Unicode value of the character at the specified location.
     * @param index The zero-based index of the desired character. If there is no character at the specified index, NaN is returned.
     */
    charCodeAt(index: number): number;

    // /**
    //  * Returns a string that contains the concatenation of two or more strings.
    //  * @param strings The strings to append to the end of the string.
    //  */
    // concat(...strings: string[]): string;

    /**
     * Returns the position of the first occurrence of a substring.
     * @param searchString The substring to search for in the string
     * @param position The index at which to begin searching the String object. If omitted, search starts at the beginning of the string.
     */
    indexOf(searchString: string, position?: number): number;

    // /**
    //  * Returns the last occurrence of a substring in the string.
    //  * @param searchString The substring to search for.
    //  * @param position The index at which to begin searching. If omitted, the search begins at the end of the string.
    //  */
    // lastIndexOf(searchString: string, position?: number): number;

    // /**
    //  * Determines whether two strings are equivalent in the current locale.
    //  * @param that String to compare to target string
    //  */
    // localeCompare(that: string): number;

    // /**
    //  * Matches a string with a regular expression, and returns an array containing the results of that search.
    //  * @param regexp A variable name or string literal containing the regular expression pattern and flags.
    //  */
    // match(regexp: string | RegExp): RegExpMatchArray | null;

    // /**
    //  * Replaces text in a string, using a regular expression or search string.
    //  * @param searchValue A string or regular expression to search for.
    //  * @param replaceValue A string containing the text to replace. When the {@linkcode searchValue} is a `RegExp`, all matches are replaced if the `g` flag is set (or only those matches at the beginning, if the `y` flag is also present). Otherwise, only the first match of {@linkcode searchValue} is replaced.
    //  */
    // replace(searchValue: string | RegExp, replaceValue: string): string;

    // /**
    //  * Replaces text in a string, using a regular expression or search string.
    //  * @param searchValue A string to search for.
    //  * @param replacer A function that returns the replacement text.
    //  */
    // replace(searchValue: string | RegExp, replacer: (substring: string, ...args: any[]) => string): string;

    // /**
    //  * Finds the first substring match in a regular expression search.
    //  * @param regexp The regular expression pattern and applicable flags.
    //  */
    // search(regexp: string | RegExp): number;

    // /**
    //  * Returns a section of a string.
    //  * @param start The index to the beginning of the specified portion of stringObj.
    //  * @param end The index to the end of the specified portion of stringObj. The substring includes the characters up to, but not including, the character indicated by end.
    //  * If this value is not specified, the substring continues to the end of stringObj.
    //  */
    // slice(start?: number, end?: number): string;

    // /**
    //  * Split a string into substrings using the specified separator and return them as an array.
    //  * @param separator A string that identifies character or characters to use in separating the string. If omitted, a single-element array containing the entire string is returned.
    //  * @param limit A value used to limit the number of elements returned in the array.
    //  */
    // split(separator: string | RegExp, limit?: number): string[];

    /**
     * Returns the substring at the specified location within a String object.
     * @param start The zero-based index number indicating the beginning of the substring.
     * @param end Zero-based index number indicating the end of the substring. The substring includes the characters up to, but not including, the character indicated by end.
     * If end is omitted, the characters from start through the end of the original string are returned.
     */
    substring(start: number, end?: number): string;

    // /** Converts all the alphabetic characters in a string to lowercase. */
    // toLowerCase(): string;

    // /** Converts all alphabetic characters to lowercase, taking into account the host environment's current locale. */
    // toLocaleLowerCase(locales?: string | string[]): string;

    // /** Converts all the alphabetic characters in a string to uppercase. */
    // toUpperCase(): string;

    // /** Returns a string where all alphabetic characters have been converted to uppercase, taking into account the host environment's current locale. */
    // toLocaleUpperCase(locales?: string | string[]): string;

    // /** Removes the leading and trailing white space and line terminator characters from a string. */
    // trim(): string;

    /** Returns the length of a String object. */
    readonly length: number;

    // // IE extensions
    // /**
    //  * Gets a substring beginning at the specified location and having the specified length.
    //  * @deprecated A legacy feature for browser compatibility
    //  * @param from The starting position of the desired substring. The index of the first character in the string is zero.
    //  * @param length The number of characters to include in the returned substring.
    //  */
    // substr(from: number, length?: number): string;

    // /** Returns the primitive value of the specified object. */
    // valueOf(): string;

    // // from typescript/lib/lib.es2015.core.d.ts

    // /**
    //  * Returns a nonnegative integer Number less than 1114112 (0x110000) that is the code point
    //  * value of the UTF-16 encoded code point starting at the string element at position pos in
    //  * the String resulting from converting this object to a String.
    //  * If there is no element at that position, the result is undefined.
    //  * If a valid UTF-16 surrogate pair does not begin at pos, the result is the code unit at pos.
    //  */
    // codePointAt(pos: number): number | undefined;

    // /**
    //  * Returns true if searchString appears as a substring of the result of converting this
    //  * object to a String, at one or more positions that are
    //  * greater than or equal to position; otherwise, returns false.
    //  * @param searchString search string
    //  * @param position If position is undefined, 0 is assumed, so as to search all of the String.
    //  */
    // includes(searchString: string, position?: number): boolean;

    // /**
    //  * Returns true if the sequence of elements of searchString converted to a String is the
    //  * same as the corresponding elements of this object (converted to a String) starting at
    //  * endPosition – length(this). Otherwise returns false.
    //  */
    // endsWith(searchString: string, endPosition?: number): boolean;

    // /**
    //  * Returns the String value result of normalizing the string into the normalization form
    //  * named by form as specified in Unicode Standard Annex #15, Unicode Normalization Forms.
    //  * @param form Applicable values: "NFC", "NFD", "NFKC", or "NFKD", If not specified default
    //  * is "NFC"
    //  */
    // normalize(form: "NFC" | "NFD" | "NFKC" | "NFKD"): string;

    // /**
    //  * Returns the String value result of normalizing the string into the normalization form
    //  * named by form as specified in Unicode Standard Annex #15, Unicode Normalization Forms.
    //  * @param form Applicable values: "NFC", "NFD", "NFKC", or "NFKD", If not specified default
    //  * is "NFC"
    //  */
    // normalize(form?: string): string;

    // /**
    //  * Returns a String value that is made from count copies appended together. If count is 0,
    //  * the empty string is returned.
    //  * @param count number of copies to append
    //  */
    // repeat(count: number): string;

    // /**
    //  * Returns true if the sequence of elements of searchString converted to a String is the
    //  * same as the corresponding elements of this object (converted to a String) starting at
    //  * position. Otherwise returns false.
    //  */
    // startsWith(searchString: string, position?: number): boolean;

    // /**
    //  * Returns an `<a>` HTML anchor element and sets the name attribute to the text value
    //  * @deprecated A legacy feature for browser compatibility
    //  * @param name
    //  */
    // anchor(name: string): string;

    // /**
    //  * Returns a `<big>` HTML element
    //  * @deprecated A legacy feature for browser compatibility
    //  */
    // big(): string;

    // /**
    //  * Returns a `<blink>` HTML element
    //  * @deprecated A legacy feature for browser compatibility
    //  */
    // blink(): string;

    // /**
    //  * Returns a `<b>` HTML element
    //  * @deprecated A legacy feature for browser compatibility
    //  */
    // bold(): string;

    // /**
    //  * Returns a `<tt>` HTML element
    //  * @deprecated A legacy feature for browser compatibility
    //  */
    // fixed(): string;

    // /**
    //  * Returns a `<font>` HTML element and sets the color attribute value
    //  * @deprecated A legacy feature for browser compatibility
    //  */
    // fontcolor(color: string): string;

    // /**
    //  * Returns a `<font>` HTML element and sets the size attribute value
    //  * @deprecated A legacy feature for browser compatibility
    //  */
    // fontsize(size: number): string;

    // /**
    //  * Returns a `<font>` HTML element and sets the size attribute value
    //  * @deprecated A legacy feature for browser compatibility
    //  */
    // fontsize(size: string): string;

    // /**
    //  * Returns an `<i>` HTML element
    //  * @deprecated A legacy feature for browser compatibility
    //  */
    // italics(): string;

    // /**
    //  * Returns an `<a>` HTML element and sets the href attribute value
    //  * @deprecated A legacy feature for browser compatibility
    //  */
    // link(url: string): string;

    // /**
    //  * Returns a `<small>` HTML element
    //  * @deprecated A legacy feature for browser compatibility
    //  */
    // small(): string;

    // /**
    //  * Returns a `<strike>` HTML element
    //  * @deprecated A legacy feature for browser compatibility
    //  */
    // strike(): string;

    // /**
    //  * Returns a `<sub>` HTML element
    //  * @deprecated A legacy feature for browser compatibility
    //  */
    // sub(): string;

    // /**
    //  * Returns a `<sup>` HTML element
    //  * @deprecated A legacy feature for browser compatibility
    //  */
    // sup(): string;

    // // from typescript/lib/lib.es2017.string.d.ts

    // /**
    //  * Pads the current string with a given string (possibly repeated) so that the resulting string reaches a given length.
    //  * The padding is applied from the start (left) of the current string.
    //  *
    //  * @param maxLength The length of the resulting string once the current string has been padded.
    //  *        If this parameter is smaller than the current string's length, the current string will be returned as it is.
    //  *
    //  * @param fillString The string to pad the current string with.
    //  *        If this string is too long, it will be truncated and the left-most part will be applied.
    //  *        The default value for this parameter is " " (U+0020).
    //  */
    // padStart(maxLength: number, fillString?: string): string;

    // /**
    //  * Pads the current string with a given string (possibly repeated) so that the resulting string reaches a given length.
    //  * The padding is applied from the end (right) of the current string.
    //  *
    //  * @param maxLength The length of the resulting string once the current string has been padded.
    //  *        If this parameter is smaller than the current string's length, the current string will be returned as it is.
    //  *
    //  * @param fillString The string to pad the current string with.
    //  *        If this string is too long, it will be truncated and the left-most part will be applied.
    //  *        The default value for this parameter is " " (U+0020).
    //  */
    // padEnd(maxLength: number, fillString?: string): string;

    // // from typescript/lib/lib.es2019.string.d.ts

    // /** Removes the trailing white space and line terminator characters from a string. */
    // trimEnd(): string;

    // /** Removes the leading white space and line terminator characters from a string. */
    // trimStart(): string;

    // /**
    //  * Removes the leading white space and line terminator characters from a string.
    //  * @deprecated A legacy feature for browser compatibility. Use `trimStart` instead
    //  */
    // trimLeft(): string;

    // /**
    //  * Removes the trailing white space and line terminator characters from a string.
    //  * @deprecated A legacy feature for browser compatibility. Use `trimEnd` instead
    //  */
    // trimRight(): string;

    // // from typescript/lib/lib.es2020.string.d.ts

    // /**
    //  * Matches a string with a regular expression, and returns an iterable of matches
    //  * containing the results of that search.
    //  * @param regexp A variable name or string literal containing the regular expression pattern and flags.
    //  */
    // matchAll(regexp: RegExp): RegExpStringIterator<RegExpExecArray>;

    // /** Converts all alphabetic characters to lowercase, taking into account the host environment's current locale. */
    // toLocaleLowerCase(locales?: Intl.LocalesArgument): string;

    // /** Returns a string where all alphabetic characters have been converted to uppercase, taking into account the host environment's current locale. */
    // toLocaleUpperCase(locales?: Intl.LocalesArgument): string;

    // /**
    //  * Determines whether two strings are equivalent in the current or specified locale.
    //  * @param that String to compare to target string
    //  * @param locales A locale string or array of locale strings that contain one or more language or locale tags. If you include more than one locale string, list them in descending order of priority so that the first entry is the preferred locale. If you omit this parameter, the default locale of the JavaScript runtime is used. This parameter must conform to BCP 47 standards; see the Intl.Collator object for details.
    //  * @param options An object that contains one or more properties that specify comparison options. see the Intl.Collator object for details.
    //  */
    // localeCompare(that: string, locales?: Intl.LocalesArgument, options?: Intl.CollatorOptions): number;
}


export type StringLike = string | String | CustomString

export function isStringLike(s: unknown): s is StringLike {
    return typeof s === 'string' || s instanceof String || (s as CustomString).__string_like__;
}

/** Try to convert `StringLike` to a primitive `string`. Might fail if the contents is longer that max allowed string length. */
export function stringLikeToString(s: StringLike): string {
    try {
        return s.toString();
    } catch (err) {
        throw new Error(`Failed to convert StringLike object into string. This might be because the length ${s.length} exceeds maximum allowed string length ${MAX_STRING_LENGTH}. (${err})`);
    }
}

/** Maximum allowed string length (might be bigger for some engines, but in Chrome and Node it is this). */
const MAX_STRING_LENGTH = 536_870_888;

/** String chunk size to be used by `ChunkedBigString`. */
const STRING_CHUNK_SIZE = 2 ** 28; // largest power of 2 smaller than MAX_STRING_LENGTH
// const STRING_CHUNK_SIZE = 8; // DEBUG, TODO use a sensible value, close to MAX_STRING_LENGTH (e.g. 2**28)


export class ChunkedBigString implements CustomString {
    readonly __string_like__: true;
    private _chunks: string[] = [];

    readonly [index: number]: string; // implemented in the constructor

    constructor() {
        return new Proxy(this, {
            // overrides getting this[index]
            get(self, index) {
                if (typeof index === 'symbol') return self[index as any];
                const i = Number.parseInt(index);
                if (isNaN(i)) return self[index as any];
                else return self.at(i);
                // TODO See how this proxy thing affects performance (CIF parsing itself uses custom parseInt because Number.parseInt is slow)
                // and possibly avoid it in favor of using .at() everywhere instead of [].
            },
        });
    }

    static fromString(content: string): ChunkedBigString {
        const out = new ChunkedBigString();
        chopString(content, STRING_CHUNK_SIZE, out._chunks);
        return out;
    }

    // potentially more performant impl, but ugly code, TODO benchmark and choose
    static fromStrings_(content: string[]): ChunkedBigString {
        const out = new ChunkedBigString();
        shakeStringChunks(content, STRING_CHUNK_SIZE, out._chunks);
        return out;
    }

    static fromStrings(content: string[]): ChunkedBigString {
        const out = new ChunkedBigString();
        for (const inputChunk of content) {
            out._append(inputChunk);
        }
        return out;
    }

    /** This is for NodeJs only, browsers don't have `Buffer`, methinks */
    static fromUtf8Buffer(buffer: Buffer, start: number = 0, end: number = buffer.length): ChunkedBigString {
        const bufferChunkSize = STRING_CHUNK_SIZE; // n bytes will always decode to <=n characters
        console.log('ChunkedBigString.fromUtf8Buffer length', buffer.length, 'isAscii', buffer.every(x => x < 128))
        const stringChunks: string[] = [];
        let readStart = start;
        while (readStart < end) {
            // TODO ensure multi-byte characters on chunk boundary are parsed correctly!
            let readEnd = Math.min(readStart + bufferChunkSize, end);
            if (readEnd < end) {
                // this is chunk boundary, adjust to avoid cutting multi-byte characters
                while ((buffer[readEnd] & 0xC0) === 0x80) { // byte after the cut is a continuation byte (10xxxxxx)
                    readEnd--;
                }
                if (readEnd === readStart) throw new Error('Input is rubbish');
            } // else this is end of read region -> default error handling
            console.log('decoding', readStart, '-', readEnd, '=', readEnd - readStart)
            const stringChunk = buffer.toString('utf-8', readStart, readEnd);
            stringChunks.push(stringChunk);
            readStart = readEnd;
        }
        return ChunkedBigString.fromStrings(stringChunks);
    }

    private _append(inputChunk: string): void {
        const chunkSize = STRING_CHUNK_SIZE;
        const tail = (this._chunks.length === 0 || this._chunks[this._chunks.length - 1].length === chunkSize) ? '' : this._chunks.pop()!;
        if (tail.length + inputChunk.length <= chunkSize) {
            this._chunks.push(tail + inputChunk);
        } else {
            let inputPtr = chunkSize - tail.length;
            this._chunks.push(tail + inputChunk.substring(0, inputPtr));
            for (; inputPtr < inputChunk.length; inputPtr += chunkSize) {
                this._chunks.push(inputChunk.substring(inputPtr, inputPtr + chunkSize));
            }
        }
        // TODO optimize - avoid concat when tail==='' (expected on ASCII inputs)
    }

    private getPosition(index: number): [iChunk: number, indexInChunk: number, outOfRange: boolean] {
        const iChunk = Math.floor(index / STRING_CHUNK_SIZE);
        const indexInChunk = index - iChunk * STRING_CHUNK_SIZE;
        const outOfRange = index < 0 || index >= this.length;
        return [iChunk, indexInChunk, outOfRange];
        // TODO think about optimizations to avoid creating array
    }

    at(index: number): string | undefined {
        const [iChunk, indexInChunk, outOfRange] = this.getPosition(index);
        if (outOfRange) return undefined;
        return this._chunks[iChunk][indexInChunk];
    }

    toString(): string {
        try {
            return this._chunks.join('');
        } catch (err) {
            throw new Error(`Failed to convert StringLike object into string. This might be because the length ${this.length} exceeds maximum allowed string length ${MAX_STRING_LENGTH}. (${err})`);
        }
    }

    _repr(): string {
        return this._chunks.join('-');
    }

    charCodeAt(index: number): number {
        const [iChunk, indexInChunk, outOfRange] = this.getPosition(index);
        if (outOfRange) return NaN;
        return this._chunks[iChunk].charCodeAt(indexInChunk);
    }

    indexOf(searchString: string, position?: number): number {
        throw new Error('NotImplementedError');
    }

    substring(start?: number, end?: number): string { // optional `start` not part of contract but works in Chrome
        const start_ = Math.min(Math.max(start ?? 0, 0), this.length);
        const end_ = Math.min(Math.max(end ?? this.length, 0), this.length);
        if (end_ < start_) {
            return this.substring(end_, start_);
        }

        if (end_ - start_ > MAX_STRING_LENGTH) {
            throw new Error(`Trying to create get a substring longer (${end_ - start_}) than maximum allowed string length (${MAX_STRING_LENGTH}).`);
        }

        const [iFirstChunk, indexInChunkFrom] = this.getPosition(start_);
        const [iLastChunk, indexInChunkTo] = this.getPosition(end_);
        const newChunks: string[] = [];
        if (iFirstChunk === iLastChunk) {
            newChunks.push(this._chunks[iFirstChunk].substring(indexInChunkFrom, indexInChunkTo));
        } else {
            newChunks.push(this._chunks[iFirstChunk].substring(indexInChunkFrom, STRING_CHUNK_SIZE));
            for (let iChunk = iFirstChunk + 1; iChunk < iLastChunk; iChunk++) {
                newChunks.push(this._chunks[iChunk]);
            }
            newChunks.push(this._chunks[iLastChunk].substring(0, indexInChunkTo));
        }
        return newChunks.join('');
    }

    get length(): number {
        const nChunks = this._chunks.length;
        if (nChunks === 0) return 0;
        return (nChunks - 1) * STRING_CHUNK_SIZE + this._chunks[nChunks - 1].length;
        // TODO compute at creation
    }
}

function chopString(content: string, chunkSize: number, out: string[]): void {
    const nChunks = Math.ceil(content.length / chunkSize);
    for (let i = 0; i < nChunks; i++) {
        out.push(content.substring(i * chunkSize, (i + 1) * chunkSize));
    }
}

function shakeStringChunks(input: string[], chunkSize: number, out: string[]): void {
    let current = '';
    for (const inputChunk of input) {
        if (current.length + inputChunk.length <= chunkSize) {
            current += inputChunk;
            if (current.length === chunkSize) {
                out.push(current);
                current = '';
            }
        } else {
            let inputPtr = chunkSize - current.length;
            current += inputChunk.substring(0, inputPtr);
            out.push(current);
            for (; inputPtr < inputChunk.length; inputPtr += chunkSize) {
                current = inputChunk.substring(inputPtr, Math.min(inputPtr + chunkSize, inputChunk.length));
                if (current.length === chunkSize) {
                    out.push(current);
                    current = '';
                }
            }
        }
    }
    if (current.length > 0) {
        out.push(current);
    }
    // TODO write proper tests
}

export function decodeBigUtf8String(buffer: Buffer, start: number = 0, end: number = buffer.length): StringLike {
    // return ChunkedBigString.fromUtf8Buffer(buffer, start, end); // DEBUG
    if (end - start <= MAX_STRING_LENGTH) {
        return buffer.toString('utf-8', start, end);
    }
    const out = ChunkedBigString.fromUtf8Buffer(buffer, start, end);
    return out.length <= MAX_STRING_LENGTH ? out.toString() : out;
}
