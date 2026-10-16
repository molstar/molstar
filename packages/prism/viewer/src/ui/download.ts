/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The one way the viewer saves a file through the page: operations return bytes (a Blob and a file name), and the
// views and commands hand them here. An object URL on an `<a download>`, released after the click.

export function downloadBlob(blob: Blob, filename: string): void {
  const url = URL.createObjectURL(blob);
  const a = document.createElement('a');
  a.href = url;
  a.download = filename;
  a.style.display = 'none';
  document.body.appendChild(a);
  a.click();
  a.remove();
  // The browser starts the download on the click; the URL is released a moment later.
  setTimeout(() => URL.revokeObjectURL(url), 1000);
}
