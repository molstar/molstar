/**
 * Copyright (c) 2024 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Adam Midlik <midlik@gmail.com>
 */


// Using dodgy typing to avoid importing peer dependency directly:
type gsStorage = any; // import { type Storage as gsStorage } from '@google-cloud/storage';
interface gsCreateReadStreamOptions {
    userProject?: string;
    validation?: 'md5' | 'crc32c' | false | true;
    start?: number;
    end?: number;
    decompress?: boolean;
} // import { type CreateReadStreamOptions as gsCreateReadStreamOptions } from '@google-cloud/storage';

let _gsClient: gsStorage | undefined;

function getGsClient(): gsStorage {
    if (!_gsClient) {
        const gsModule = getGsModule();
        _gsClient = new gsModule.Storage();
    }
    return _gsClient;
}

function getGsModule() {
    try {
        return require('@google-cloud/storage');
    } catch (err) {
        console.error('Dynamic import of "@google-cloud/storage" failed. If you want to use ModelServer/VolumeServer with Google Cloud Storage, install "@google-cloud/storage" package (see peerDependencies in package.json).');
        throw err;
    }
}


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

export async function downloadGs(gsUrl: string, options?: gsCreateReadStreamOptions): Promise<Buffer> {
    const { bucket, file } = parseGsUrl(gsUrl);
    try {
        const response = await getGsClient().bucket(bucket).file(file).download(options);
        const buffer = response[0];
        return buffer;
    } catch (err) {
        if (err.code === 401) {
            console.error('Error 401: Unauthorized access to Google Cloud Storage. To set up authorization, run `export GOOGLE_APPLICATION_CREDENTIALS=/path/to/credential/file.json`');
        }
        throw err;
    }
}
