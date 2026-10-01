/**
 * Copyright (c) 2025 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

export async function fileToDataUri(file: File): Promise<string> {
    const filename = file.name.toLowerCase() || 'file';
    const isImage = ['jpg', 'jpeg', 'png', 'gif', 'webp'].some(ext => filename.endsWith(`.${ext}`));
    const isAudio = ['mp3', 'wav', 'ogg'].some(ext => filename.endsWith(`.${ext}`));

    let type = 'application/octet-stream';
    if (isImage) {
        const ext = filename.split('.').pop()?.toLowerCase();
        switch (ext) {
            case 'jpg':
                type = 'image/jpeg';
                break;
            default:
                type = `image/${ext}`;
                break;
        }
    } else if (isAudio) {
        const ext = filename.split('.').pop()?.toLowerCase();
        switch (ext) {
            case 'mp3':
                type = 'audio/mpeg';
                break;
            default:
                type = `audio/${ext}`;
                break;
        }
    }

    const bytes = await file.arrayBuffer();
    const reader = new FileReader();
    reader.readAsDataURL(new Blob([bytes], { type }));
    const data = await new Promise<string>((resolve, reject) => {
        reader.onload = () => resolve(reader.result as string);
        reader.onerror = () => reject(reader.error);
    });

    return data;
}