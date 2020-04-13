/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

export function readFile (file: File, isBinary = false) {
    const fileReader = new FileReader();
    return new Promise<string | Uint8Array>((resolve, reject) => {
        fileReader.onerror = () => {
            fileReader.abort();
            reject(new DOMException('Error parsing file.'));
        };
        fileReader.onload = () => {
            resolve(isBinary ? new Uint8Array(fileReader.result as ArrayBuffer) : fileReader.result as string);
        };
        if (isBinary) {
            fileReader.readAsArrayBuffer(file);
        } else {
            fileReader.readAsText(file);
        }
    });
}

export function readFileAsText(file: File) {
    return readFile(file, false) as Promise<string>;
}

export function readFileAsBuffer(file: File) {
    return readFile(file, true) as Promise<Uint8Array>;
}

export async function readUrl(url: string, isBinary: boolean) {
    const response = await fetch(url);
    return isBinary ? new Uint8Array(await response.arrayBuffer()) : await response.text();
}

export function readUrlAsText(url: string) {
    return readUrl(url, false) as Promise<string>;
}

export function readUrlAsBuffer(url: string) {
    return readUrl(url, true) as Promise<Uint8Array>;
}