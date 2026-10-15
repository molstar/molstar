/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

/** Only lossless for opaque images: 2D canvases store premultiplied alpha */
export function imageDataToCanvas(image: ImageData) {
  const canvas = document.createElement('canvas');
  canvas.width = image.width;
  canvas.height = image.height;
  canvas.getContext('2d')!.putImageData(image, 0, 0);
  return canvas;
}

export function encodePng(image: ImageData): Promise<Blob> {
  return new Promise((resolve, reject) => {
    imageDataToCanvas(image).toBlob((b) => (b ? resolve(b) : reject(new Error('PNG encoding failed'))), 'image/png');
  });
}

export async function encodePngDataUrl(image: ImageData): Promise<string> {
  const blob = await encodePng(image);
  return new Promise((resolve, reject) => {
    const reader = new FileReader();
    reader.onload = () => resolve(reader.result as string);
    reader.onerror = () => reject(reader.error);
    reader.readAsDataURL(blob);
  });
}

/** Returns undefined if the image does not exist */
export async function loadPng(url: string): Promise<ImageData | undefined> {
  const response = await fetch(url, { cache: 'no-store' });
  if (!response.ok) return undefined;
  const bitmap = await createImageBitmap(await response.blob(), {
    premultiplyAlpha: 'none',
    colorSpaceConversion: 'none',
  });
  const canvas = document.createElement('canvas');
  canvas.width = bitmap.width;
  canvas.height = bitmap.height;
  const ctx = canvas.getContext('2d', { willReadFrequently: true })!;
  ctx.drawImage(bitmap, 0, 0);
  bitmap.close();
  return ctx.getImageData(0, 0, canvas.width, canvas.height);
}

export async function downloadPng(image: ImageData, fileName: string) {
  const url = URL.createObjectURL(await encodePng(image));
  const a = document.createElement('a');
  a.href = url;
  a.download = fileName;
  a.click();
  setTimeout(() => URL.revokeObjectURL(url), 1000);
}
