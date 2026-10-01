/*
 * Copyright (c) 2017 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * from https://github.com/dsehnal/CIFTools.js
 * @author David Sehnal <david.sehnal@gmail.com>
 */

type ReaderResult<T> = ReaderResult.Success<T> | ReaderResult.Error

namespace ReaderResult {
    export function error<T>(message: string, line = -1): ReaderResult<T> {
        return new Error(message, line);
    }

    export function success<T>(result: T, warnings: string[] = []): ReaderResult<T> {
        return new Success<T>(result, warnings);
    }

    export class Error {
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

    export class Success<T> {
        isError: false = false;

        constructor(public result: T, public warnings: string[]) { }
    }
}

export { ReaderResult };