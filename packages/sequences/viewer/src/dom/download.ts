/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The UI's download (seqspec:interaction#the-default-ui, `Export`): `offer` for `browserExportAction` in the registry
// of a viewer the element makes.

/** The finished SVG as a file download: a Blob URL on a temporary link, released a little later. */
export function downloadSvg(svg: string, name: string): void {
  const url = URL.createObjectURL(new Blob([svg], { type: 'image/svg+xml' }));
  const link = Object.assign(document.createElement('a'), { href: url, download: name });
  link.click();
  setTimeout(() => URL.revokeObjectURL(url), 10_000);
}
