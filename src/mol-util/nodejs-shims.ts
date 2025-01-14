/**
 * Copyright (c) 2018-2023 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Adam Midlik <midlik@gmail.com>
 * @author Russell Parker <russell@benchling.com>
 *
 * Implements some browser-only global variables for Node.js environment.
 * These workarounds will also work in browsers as usual.
 */


/** Determines whether the current code is running in Node.js */
export const RUNNING_IN_NODEJS = typeof process !== 'undefined' && process.versions != null && process.versions.node != null;

/** Like `XMLHttpRequest` but works also in Node.js */
export const XMLHttpRequest_ = getXMLHttpRequest();

/** Like `File` but works also in Node.js */
export const File_ = getFile();


function getXMLHttpRequest(): typeof XMLHttpRequest {
    if (typeof XMLHttpRequest === 'undefined' || RUNNING_IN_NODEJS) {
        return require('xhr2');
    } else {
        return XMLHttpRequest;
    }
}

function getFile(): typeof File {
    if (typeof File === 'undefined' || RUNNING_IN_NODEJS) {
        class File_NodeJs implements File {
            private readonly blob: Blob;
            // Blob fields
            readonly size: number;
            readonly type: string;
            arrayBuffer() { return this.blob.arrayBuffer(); }
            slice(start?: number, end?: number, contentType?: string) { return this.blob.slice(start, end, contentType); }
            stream() { return this.blob.stream(); }
            text() { return this.blob.text(); }
            bytes() { return this.blob.bytes(); }
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
