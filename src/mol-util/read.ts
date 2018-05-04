/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

export const readFileAs = (file: File, isBinary = false) => {
    const fileReader = new FileReader()
    return new Promise<string | Uint8Array>((resolve, reject) => {
        fileReader.onerror = () => {
            fileReader.abort()
            reject(new DOMException('Error parsing file.'))
        }
        fileReader.onload = () => {
            resolve(isBinary ? new Uint8Array(fileReader.result) : fileReader.result)
        }
        if (isBinary) {
            fileReader.readAsArrayBuffer(file)
        } else {
            fileReader.readAsText(file)
        }
    })
}

export function readFileAsText(file: File) {
    return readFileAs(file, false) as Promise<string>
}

export function readFileAsBuffer(file: File) {
    return readFileAs(file, true) as Promise<Uint8Array>
}

export async function readUrlAs(url: string, isBinary: boolean) {
    const response = await fetch(url);
    return isBinary ? new Uint8Array(await response.arrayBuffer()) : await response.text();
}

export function readUrlAsText(url: string) {
    return readUrlAs(url, false) as Promise<string>
}

export function readUrlAsBuffer(url: string) {
    return readUrlAs(url, true) as Promise<Uint8Array>
}