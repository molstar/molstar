/**
 * Copyright (c) 2024 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Adam Midlik <midlik@gmail.com>
 */

import { DownloadOptions, Storage } from '@google-cloud/storage';


function parseGsUrl(url: string) {
    if (!url.startsWith('gs://')) throw new Error(`Invalid Google Cloud Storage URL (must start with 'gs://'): ${url}`);
    const bucket_file = url.slice('gs://'.length);
    const slashIndex = bucket_file.indexOf('/');
    if (slashIndex < 0) throw new Error(`Invalid Google Cloud Storage URL (missing slash): ${url}`);
    return {
        bucket: bucket_file.slice(0, slashIndex),
        file: bucket_file.slice(slashIndex + 1, undefined),
    };
}

let _gsClient: Storage | undefined;
function getGsClient() {
    return _gsClient ??= new Storage();
}

export async function downloadGs(gsUrl: string, options?: DownloadOptions) {
    const { bucket, file } = parseGsUrl(gsUrl);
    const response = await getGsClient().bucket(bucket).file(file).download(options);
    const buffer = response[0];
    return buffer;
}
