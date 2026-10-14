/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The SVG as a PNG, by Chrome (`playwright`, the Chrome the browser tests use): the page shows the SVG at its own size and is
// photographed. Node has no rasterizer of its own; only `seqv-render --png` and tests that want pixels import this, and
// a CLI user without playwright or Chrome is told which is missing.

const first = (e: unknown) => String(e instanceof Error ? e.message : e).split('\n')[0];

/** Chrome through playwright, or an error that names what is missing. */
async function chrome() {
  let pw: typeof import('playwright');
  try {
    pw = await import('playwright');
  } catch (e) {
    throw new Error(`--png needs the \`playwright\` package (npm install playwright) and Google Chrome: ${first(e)}`);
  }
  try {
    return await pw.chromium.launch({ channel: 'chrome' });
  } catch (e) {
    throw new Error(`--png needs Google Chrome (npx playwright install chrome): ${first(e)}`);
  }
}

/** `svg` rasterized at the size its `width` and `height` give (its scale is in them), on a transparent background. */
export async function rasterizeSvg(svg: string): Promise<Uint8Array> {
  const size = /<svg[^>]* width="([\d.]+)" height="([\d.]+)"/.exec(svg);
  if (!size) throw new Error('rasterizeSvg: the text is not an SVG with a width and a height.');
  const browser = await chrome();
  try {
    const page = await browser.newPage({
      viewport: { width: Math.ceil(Number(size[1])), height: Math.ceil(Number(size[2])) },
    });
    await page.setContent(`<!doctype html><body style="margin:0">${svg}</body>`);
    return await page.screenshot({ omitBackground: true });
  } finally {
    await browser.close();
  }
}
