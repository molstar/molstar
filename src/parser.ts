/*
 * Copyright (c) 2017 molio contributors, licensed under MIT, See LICENSE file for more info.
 *
 * from https://github.com/dsehnal/CIFTools.js
 * @author David Sehnal <david.sehnal@gmail.com>
 */

export type ParserResult<T> = ParserSuccess<T> | ParserError

export namespace ParserResult {
    export function error<T>(message: string, line = -1): ParserResult<T> {
        return new ParserError(message, line);
    }

    export function success<T>(result: T, warnings: string[] = []): ParserResult<T> {
        return new ParserSuccess<T>(result, warnings);
    }
}

export class ParserError {
    isError: true = true;

    toString() {
        if (this.line >= 0) {
            return `[Line ${this.line}] ${this.message}`;
        }
        return this.message;
    }

    constructor(
        public message: string,
        public line: number) {
    }
}

export class ParserSuccess<T> {
    isError: false = false;

    constructor(public result: T, public warnings: string[]) { }
}