/**
 * Copyright (c) 2018-2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Adam Midlik <midlik@gmail.com>
 *
 * Implements some browser-only global variables for Node.js environment.
 * These workarounds will also work in browsers as usual.
 */


/** Determines whether the current code is running in Node.js or in a browser */
export const RUNNING_IN_NODEJS = typeof document === 'undefined';

/** Like `XMLHttpRequest` but works also in Node.js */
export const XMLHttpRequest_ = getXMLHttpRequest();

/** Like `File` but works also in Node.js */
export const File_ = getFile();


function getXMLHttpRequest(): typeof XMLHttpRequest {
    if (typeof document === 'undefined') {
        return require('xhr2');
    } else {
        return XMLHttpRequest;
    }
}

function getFile(): typeof File {
    if (typeof document === 'undefined') {
        class File_NodeJs implements File {
            private readonly blob: Blob;
            // Blob fields
            readonly size: number;
            readonly type: string;
            arrayBuffer() { return this.blob.arrayBuffer(); }
            slice(start?: number, end?: number, contentType?: string) { return this.blob.slice(start, end, contentType); }
            stream() { return this.blob.stream(); }
            text() { return this.blob.text(); }
            // File fields
            name: string;
            lastModified: number;
            webkitRelativePath: string;

            constructor(fileBits: BlobPart[], fileName: string, options?: FilePropertyBag) {
                this.blob = new Blob(fileBits, options);
                // Blob fields
                this.size = this.blob.size;
                this.type = this.blob.type;
                // File fields
                this.name = fileName;
                this.lastModified = options?.lastModified ?? 0;
                this.webkitRelativePath = '';
            }
        }
        return File_NodeJs;
    } else {
        return File;
    }
}
