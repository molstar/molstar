/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import assert from 'node:assert/strict';
import { build } from 'esbuild';
import { chromium } from 'playwright';
import { createServer } from 'node:http';
import { access, mkdir, readFile, writeFile } from 'node:fs/promises';
import { resolve } from 'node:path';
import { fileURLToPath } from 'node:url';

const root = fileURLToPath(new URL('../../', import.meta.url));
const output = resolve(root, 'build/smoke-video/video.js');
await mkdir(resolve(root, 'build/smoke-video'), { recursive: true });
const bundle = await build({
  absWorkingDir: root,
  entryPoints: ['packages/extensions/video-export/src/_test/browser.ts'],
  outfile: output,
  bundle: true,
  format: 'iife',
  globalName: 'VideoSmoke',
  conditions: ['molstar-src'],
  metafile: true,
  external: ['fs', 'path', 'crypto', 'stream'],
  loader: { '.scss': 'empty', '.html': 'empty', '.sdf': 'empty' },
  define: {
    'process.env.NODE_ENV': '"production"',
    'process.env.DEBUG': 'false',
    __MOLSTAR_PLUGIN_VERSION__: '"smoke"',
  },
});
assert(
  !Object.keys(bundle.metafile.inputs).some((path) => path.includes('h264-mp4-encoder')),
  'Default browser bundle imports WASM',
);
const script = await readFile(output);
const server = createServer((req, res) => {
  res.setHeader('Content-Type', req.url === '/video.js' ? 'text/javascript' : 'text/html');
  res.end(req.url === '/video.js' ? script : '<script src="/video.js"></script>');
});
await new Promise((done, reject) => {
  server.once('error', reject);
  server.listen(0, '127.0.0.1', done);
});
let browser;
try {
  const chrome = '/Applications/Google Chrome.app/Contents/MacOS/Google Chrome';
  const executablePath =
    process.env.MOLSTAR_SMOKE_BROWSER ??
    (await access(chrome).then(
      () => chrome,
      () => undefined,
    ));
  browser = await chromium.launch({ headless: true, executablePath });
  const page = await browser.newPage();
  await page.goto(`http://127.0.0.1:${server.address().port}`);
  const results = await page.evaluate(() => window.VideoSmoke.run());
  assert(results.length > 0);
  await writeFile(resolve(root, 'build/smoke-video/results.json'), JSON.stringify(results, null, 2));
  console.log(`Mediabunny browser export passed: ${results.join(', ')}; WASM excluded from default bundle`);
} finally {
  await browser?.close();
  await new Promise((done) => server.close(done));
}
