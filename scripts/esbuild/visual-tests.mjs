import * as esbuild from 'esbuild';
import fs from 'node:fs';
import http from 'node:http';
import path from 'node:path';
import { fileURLToPath } from 'node:url';

// Bundles tests/src/<suite> into tests/build/<suite> and serves it from the repository root; the suite is `visual` (the
// visual tests, default) or `benchmark` (`--suite benchmark`).
const root = path.resolve(path.dirname(fileURLToPath(import.meta.url)), '../..');
const args = process.argv.slice(2).filter((arg) => arg !== '--');
const portIndex = args.indexOf('--port');
const port = portIndex >= 0 ? Number(args[portIndex + 1]) : 1338;
const production = args.includes('--prd');
const watch = !production && !args.includes('--no-watch');
const serve = watch || args.includes('--serve');
const suiteIndex = args.indexOf('--suite');
const suite = suiteIndex >= 0 ? args[suiteIndex + 1] : 'visual';
if (!['visual', 'benchmark'].includes(suite)) throw new Error(`unknown suite '${suite}'`);

const outdir = path.join(root, `tests/build/${suite}`);

const htmlAssetPlugin = {
  name: 'molstar-visual-tests-html',
  setup(build) {
    build.onLoad({ filter: /\.html$/ }, async ({ path: input }) => {
      await fs.promises.mkdir(outdir, { recursive: true });
      await fs.promises.copyFile(input, path.join(outdir, path.basename(input)));
      return { contents: '', loader: 'empty' };
    });
  },
};

const sourceJsExtensionPlugin = {
  name: 'molstar-source-js-extension',
  setup(build) {
    build.onResolve({ filter: /^\./ }, (resolveArgs) => {
      if (!resolveArgs.path.endsWith('.js')) return null;
      const file = path.resolve(resolveArgs.resolveDir, resolveArgs.path);
      if (fs.existsSync(file)) return null;
      for (const ext of ['.ts', '.tsx']) {
        const source = file.slice(0, -'.js'.length) + ext;
        if (fs.existsSync(source)) return { path: source };
      }
      return null;
    });
  },
};

const options = {
  absWorkingDir: root,
  entryPoints: [path.join(root, `tests/src/${suite}/index.ts`)],
  tsconfig: path.join(root, 'tests/tsconfig.json'),
  bundle: true,
  platform: 'browser',
  format: 'iife',
  outfile: path.join(outdir, 'index.js'),
  minify: production,
  minifyIdentifiers: false,
  sourcemap: !production,
  conditions: ['molstar-src', 'import', 'default'],
  external: ['crypto', 'fs', 'path', 'stream'],
  plugins: [sourceJsExtensionPlugin, htmlAssetPlugin],
  define: {
    'process.env.NODE_ENV': JSON.stringify(production ? 'production' : 'development'),
    'process.env.DEBUG': JSON.stringify(process.env.DEBUG || false),
    __MOLSTAR_PLUGIN_VERSION__: JSON.stringify(`(${suite} tests)`),
    __MOLSTAR_BUILD_TIMESTAMP__: String(Date.now()),
  },
  logLevel: 'info',
};

if (watch) {
  const ctx = await esbuild.context(options);
  await ctx.watch();
} else {
  await esbuild.build(options);
}
console.log(`${watch ? 'Watching' : 'Built'} ${suite} tests → ${path.relative(root, outdir)}`);

if (serve) {
  const mime = { '.html': 'text/html', '.js': 'text/javascript', '.map': 'application/json' };
  const server = http.createServer(async (req, res) => {
    try {
      const pathname = decodeURIComponent(new URL(req.url, 'http://localhost').pathname);
      const resolved = path.resolve(root, `.${pathname}`);
      if (!resolved.startsWith(root + path.sep)) {
        res.writeHead(403).end();
        return;
      }
      const stat = await fs.promises.stat(resolved);
      const file = stat.isDirectory() ? path.join(resolved, 'index.html') : resolved;
      const body = await fs.promises.readFile(file);
      res.writeHead(200, { 'content-type': mime[path.extname(file)] ?? 'application/octet-stream' }).end(body);
    } catch {
      res.writeHead(404).end('Not found');
    }
  });
  await new Promise((resolve) => server.listen(port, '0.0.0.0', resolve));
  console.log(`Visual tests: http://localhost:${port}/tests/build/visual/`);
}
