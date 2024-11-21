/**
 * Copyright (c) 2024 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Adam Midlik <midlik@gmail.com>
 */

import { DownloadOptions, Storage } from '@google-cloud/storage';
import fetch, { Response } from 'node-fetch';


export async function extendedFetch(url: string, init?: fetch.RequestInit): Promise<Response> {
    if (url.startsWith('gs://')) return fetch_GS(url, init);
    else return fetch(url, init);
}

async function fetch_GS(url: string, init?: fetch.RequestInit): Promise<Response> {
    const fields = parseGsUrl(url);
    const data = await downloadGs(fields.bucket, fields.file);
    return new Response(data, init);
}

export function parseGsUrl(url: string) {
    if (!url.startsWith('gs://')) throw new Error(`Invalid GS URL (must start with 'gs://'): ${url}`);
    const bucket_file = url.slice('gs://'.length);
    const slashIndex = bucket_file.indexOf('/');
    if (slashIndex < 0) throw new Error(`Invalid GS URL (missing slash): ${url}`);
    return {
        bucket: bucket_file.slice(0, slashIndex),
        file: bucket_file.slice(slashIndex + 1, undefined),
    };
}

let _gsClient: Storage | undefined;
function getGsClient() {
    return _gsClient ??= new Storage();
}

export async function downloadGs(bucketName: string, srcFileName: string, options?: DownloadOptions) {
    const response = await getGsClient().bucket(bucketName).file(srcFileName).download(options); // decompress:true?, end inclusive!
    const buffer = response[0];
    return buffer;
}
