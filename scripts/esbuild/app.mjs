import * as esbuild from 'esbuild';
import fs from 'node:fs';
import path from 'node:path';
import http from 'node:http';
import https from 'node:https';
import { fileURLToPath, pathToFileURL } from 'node:url';
import { sassPlugin } from 'esbuild-sass-plugin';

const root = path.resolve(path.dirname(fileURLToPath(import.meta.url)), '../..');
const args = process.argv.slice(2).filter((arg) => arg !== '--');
const value = (flag, fallback) => {
  const i = args.indexOf(flag);
  return i >= 0 ? args[i + 1] : fallback;
};
const all = args.includes('--all');
const dev = args.includes('--dev') || args.includes('--watch');
const production = args.includes('--prd') || !dev;
// `--metafile`: esbuild's metafile beside each output (`<name>.meta.json`), for size breakdowns.
const metafile = args.includes('--metafile');
const serve = dev || args.includes('--serve');
const port = Number(value('--port', '1338'));
const kind = value('--kind');
if (kind && !['app', 'example'].includes(kind)) throw new Error('--kind must be app or example.');
const positional = [];
for (let i = 0; i < args.length; i++) {
  if (['--kind', '--port'].includes(args[i])) {
    i++;
    continue;
  }
  if (!args[i].startsWith('-')) positional.push(args[i]);
}
const version = JSON.parse(fs.readFileSync(path.join(root, 'version.json'), 'utf8')).version;
const inventoryPath = path.join(root, 'scripts/workspace/inventory.json');
const inventory = fs.existsSync(inventoryPath) ? JSON.parse(fs.readFileSync(inventoryPath, 'utf8')) : { packages: [] };
const packagesByName = new Map((inventory.packages ?? []).map((p) => [p.name, p]));
const sassImporter = {
  canonicalize(specifier, { containingUrl }) {
    let file;
    if (specifier.startsWith('pkg:')) {
      const request = specifier.slice(4);
      const packageName = request.startsWith('@')
        ? request.indexOf('/', request.indexOf('/') + 1)
        : request.indexOf('/');
      if (packageName < 0) return null;
      const pkg = packagesByName.get(request.slice(0, packageName));
      const subpath = request.slice(packageName + 1);
      if (!pkg || !subpath) return null;
      const sourceRoot = path.resolve(root, pkg.path, 'src');
      file = path.resolve(sourceRoot, subpath);
      if (!file.startsWith(sourceRoot + path.sep)) return null;
    } else if (containingUrl?.protocol === 'file:' && !specifier.startsWith('http:')) {
      file = path.resolve(path.dirname(fileURLToPath(containingUrl)), specifier);
    } else return null;
    const ext = path.extname(file);
    const candidates = ext
      ? [file]
      : [
          file + '.scss',
          path.join(path.dirname(file), `_${path.basename(file)}.scss`),
          file + '.sass',
          path.join(path.dirname(file), `_${path.basename(file)}.sass`),
        ];
    const resolved = candidates.find((candidate) => fs.existsSync(candidate));
    return resolved ? pathToFileURL(resolved) : null;
  },
  load(url) {
    const file = fileURLToPath(url);
    let contents = fs.readFileSync(file, 'utf8');
    contents = contents.replace(
      /(@(?:use|forward|import)\s+)(['"])(\.{1,2}\/[^'"]+)\2/gu,
      (match, prefix, quote, specifier) => {
        const base = path.resolve(path.dirname(file), specifier);
        const ext = path.extname(base);
        const candidates = ext
          ? [base, path.join(path.dirname(base), `_${path.basename(base)}`)]
          : [
              base + '.scss',
              path.join(path.dirname(base), `_${path.basename(base)}.scss`),
              base + '.sass',
              path.join(path.dirname(base), `_${path.basename(base)}.sass`),
            ];
        const resolved = candidates.find((candidate) => fs.existsSync(candidate));
        return resolved ? `${prefix}${quote}${pathToFileURL(resolved).href}${quote}` : match;
      },
    );
    return { contents, syntax: file.endsWith('.sass') ? 'indented' : 'scss' };
  },
};
const pkgInfo = (p) => {
  const manifest = JSON.parse(fs.readFileSync(path.join(root, p.path, 'package.json'), 'utf8'));
  return { ...p, ...manifest, config: manifest.molstar ?? {} };
};
const apps = (inventory.packages ?? []).filter((p) => {
  if (!['app', 'example'].includes(p.kind) || (kind && p.kind !== kind)) return false;
  const manifestPath = path.join(root, p.path, 'package.json');
  if (!fs.existsSync(manifestPath)) return false;
  const manifest = JSON.parse(fs.readFileSync(manifestPath, 'utf8'));
  // A package with its own builder (`molstar.builder`, a Vite example) is built by its own scripts.
  return manifest.molstar?.platform !== 'node' && !manifest.molstar?.builder;
});
const matches = (pkg, target) =>
  pkg.name === target || pkg.path === target || (kind && path.basename(pkg.path) === target);
const unknownTargets = positional.filter((target) => !apps.some((pkg) => matches(pkg, target)));
const selected =
  all || (kind && !positional.length) ? apps : apps.filter((pkg) => positional.some((target) => matches(pkg, target)));
if (args.includes('--help') || unknownTargets.length || !selected.length) {
  if (!args.includes('--help') && unknownTargets.length) console.error(`Unknown targets: ${unknownTargets.join(', ')}`);
  const usage = kind
    ? `Usage: pnpm dev:${kind}s [-- <name> ...] [--port <port>] (no names watches all)`
    : 'Usage: node scripts/esbuild/app.mjs <package-name|path> [--prd|--dev] | --all [--prd|--dev]';
  console.log(usage);
  console.log(
    `Available browser ${kind ? `${kind}s` : 'apps/examples'}: ${apps.map((p) => (kind ? path.basename(p.path) : p.path)).join(', ')}`,
  );
  process.exit(args.includes('--help') ? 0 : 2);
}
for (const p of selected)
  if (!fs.existsSync(path.join(root, p.path, 'package.json')))
    throw new Error(`Missing manifest: ${p.path}/package.json`);

function staticAssetsPlugin(outdir) {
  return {
    name: 'molstar-static-assets',
    setup(build) {
      build.onLoad({ filter: /\.(jpg|jpeg|png|gif|ico|html|htm|svg|wasm|bin|dat|sdf)$/i }, async ({ path: input }) => {
        const ext = path.extname(input).toLowerCase();
        const name = path.basename(input);
        const isImage = /\.(jpg|jpeg|png|gif|svg)$/i.test(ext);
        const destDir = isImage ? path.join(outdir, 'images') : outdir;
        await fs.promises.mkdir(destDir, { recursive: true });
        await fs.promises.copyFile(input, path.join(destDir, name));
        if (/\.(html|htm|ico|sdf)$/i.test(ext)) return { contents: '', loader: 'empty' };
        return { contents: `${isImage ? 'images/' : ''}${name}`, loader: 'text' };
      });
    },
  };
}
/** The single file's images: copied to `images/` beside it, imported as URLs relative to the module. */
function moduleRelativeImagesPlugin(outdir) {
  return {
    name: 'molstar-module-relative-images',
    setup(build) {
      build.onLoad({ filter: /\.(jpg|jpeg|png|gif|svg)$/i }, async ({ path: input }) => {
        const name = path.basename(input);
        await fs.promises.mkdir(path.join(outdir, 'images'), { recursive: true });
        await fs.promises.copyFile(input, path.join(outdir, 'images', name));
        return {
          contents: `export default new URL(${JSON.stringify(`./images/${name}`)}, import.meta.url).href;`,
          loader: 'js',
        };
      });
    },
  };
}
const sourceJsExtensionPlugin = {
  name: 'molstar-source-js-extension',
  setup(build) {
    build.onResolve({ filter: /^\./ }, (args) => {
      if (!args.path.endsWith('.js')) return null;
      const file = path.resolve(args.resolveDir, args.path);
      if (fs.existsSync(file)) return null;
      for (const ext of ['.ts', '.tsx']) {
        const source = file.slice(0, -'.js'.length) + ext;
        if (fs.existsSync(source)) return { path: source };
      }
      return null;
    });
  },
};
/**
 * zod's `locales` namespace (62 translations, ~260 kB) is re-exported whole and the AI SDKs import zod as a namespace, so
 * the assistant's chunk would carry all of them. Nothing here picks a locale: the namespace is English only (the
 * default locale zod registers itself imports `en.js` directly).
 */
const zodEnglishOnlyPlugin = {
  name: 'molstar-zod-english-only',
  setup(build) {
    build.onLoad({ filter: /[\\/]zod[\\/]v4[\\/]locales[\\/]index\.js$/ }, ({ path: input }) => ({
      contents: 'export { default as en } from "./en.js";',
      resolveDir: path.dirname(input),
      loader: 'js',
    }));
  },
};
function entryFor(pkg, config) {
  const raw = path.resolve(root, pkg.path, config.entry ?? 'src/index.ts');
  if (fs.existsSync(raw)) return raw;
  if (fs.existsSync(raw + 'x')) return raw + 'x';
  throw new Error(`${pkg.name}: configured entry does not exist: ${path.relative(root, raw)}`);
}
async function buildPackage(rawPkg) {
  const pkg = pkgInfo(rawPkg),
    cfg = pkg.config;
  const dir = path.join(root, pkg.path),
    outdir = path.join(dir, 'build');
  const entry = entryFor(pkg, cfg);
  await fs.promises.mkdir(outdir, { recursive: true });
  const filename = cfg.filename ?? (pkg.kind === 'example' ? 'index.js' : 'molstar.js');
  const outfile = path.join(outdir, filename);
  // `format: "esm"` (the Prism viewer, app.md §2): an ES module entry with lazy chunks beside it, loaded by
  // `<script type="module">`; the default is one IIFE file with a global.
  const esm = cfg.format === 'esm';
  if (esm && cfg.themes?.length) throw new Error(`${pkg.name}: themes are built as IIFE files; an esm app takes none`);
  const output = esm
    ? {
        format: 'esm',
        splitting: true,
        outdir,
        entryNames: path.basename(filename, '.js'),
        chunkNames: 'chunks/[name]-[hash]',
      }
    : { format: 'iife', globalName: cfg.globalName ?? 'molstar', outfile };
  if (esm) await fs.promises.rm(path.join(outdir, 'chunks'), { recursive: true, force: true });
  const options = {
    absWorkingDir: root,
    entryPoints: [entry],
    bundle: true,
    platform: 'browser',
    ...output,
    tsconfig: fs.existsSync(path.join(dir, 'tsconfig.json'))
      ? path.join(dir, 'tsconfig.json')
      : path.join(root, 'tsconfig.json'),
    minify: production,
    sourcemap: !production,
    conditions: ['molstar-src', 'import', 'default'],
    external: ['crypto', 'fs', 'path', 'stream'],
    plugins: [
      sourceJsExtensionPlugin,
      zodEnglishOnlyPlugin,
      staticAssetsPlugin(outdir),
      sassPlugin({ type: 'css', embedded: false, importers: [sassImporter], silenceDeprecations: ['import'] }),
    ],
    define: {
      'process.env.NODE_ENV': JSON.stringify(production ? 'production' : 'development'),
      'process.env.DEBUG': JSON.stringify(process.env.DEBUG || false),
      __MOLSTAR_PLUGIN_VERSION__: JSON.stringify(version),
    },
    logLevel: 'info',
  };
  if (pkg.kind === 'example') {
    options.plugins.push({
      name: 'molstar-example-css',
      setup(build) {
        build.onEnd(async (result) => {
          const indexCss = path.join(outdir, 'index.css');
          if (!result.errors.length && fs.existsSync(indexCss)) {
            await fs.promises.rename(indexCss, path.join(outdir, 'molstar.css'));
          }
        });
      },
    });
  }
  if (dev) {
    const ctx = await esbuild.context(options);
    await ctx.watch();
    contexts.push(ctx);
  } else {
    const result = await esbuild.build({ ...options, metafile });
    if (metafile)
      await fs.promises.writeFile(
        path.join(outdir, `${path.basename(filename, '.js')}.meta.json`),
        JSON.stringify(result.metafile),
      );
  }
  // `single` (the single-file Prism viewer, app.md §2): one more entry built as one ES module, its lazy parts inlined
  // (no splitting: esbuild bundles dynamic imports into the file), so a page imports it with no other request.
  if (cfg.single) {
    const singleOptions = {
      ...options,
      entryPoints: [path.resolve(dir, cfg.single.entry)],
      format: 'esm',
      splitting: false,
      outdir: undefined,
      entryNames: undefined,
      chunkNames: undefined,
      outfile: path.join(outdir, cfg.single.filename),
      // Images beside the file (`images/`), by a URL relative to the module (not the page): fetched only when used (a
      // background's picture), instead of inlined into every page's download.
      plugins: [
        ...options.plugins.filter((p) => p.name !== 'molstar-static-assets'),
        moduleRelativeImagesPlugin(outdir),
      ],
    };
    if (dev) {
      const ctx = await esbuild.context(singleOptions);
      await ctx.watch();
      contexts.push(ctx);
    } else {
      const result = await esbuild.build({ ...singleOptions, metafile });
      if (metafile) {
        const name = path.basename(cfg.single.filename, '.js');
        await fs.promises.writeFile(path.join(outdir, `${name}.meta.json`), JSON.stringify(result.metafile));
      }
    }
  }
  // Size budgets (app.md §1): a production build fails when a file outgrows its budget (bytes, `molstar.budgets`).
  if (!dev) {
    for (const f of fs.readdirSync(outdir)) if (f.endsWith('.js.map')) await fs.promises.rm(path.join(outdir, f));
    const over = Object.entries(cfg.budgets ?? {})
      .map(([file, budget]) => [file, budget, fs.statSync(path.join(outdir, file)).size])
      .filter(([, budget, size]) => size > budget);
    if (over.length) {
      throw new Error(
        `${pkg.name}: over budget: ${over.map(([f, b, s]) => `${f} ${(s / 1e6).toFixed(2)} MB > ${(b / 1e6).toFixed(2)} MB`).join(', ')}`,
      );
    }
  }
  for (const theme of cfg.themes ?? []) {
    const themeRoot = path.join(dir, 'src/theme');
    let themeEntry = path.join(themeRoot, `${theme}.ts`);
    if (!fs.existsSync(themeEntry)) themeEntry += 'x';
    if (!fs.existsSync(themeEntry)) throw new Error(`${pkg.name}: theme entry not found for ${theme}`);
    const themeOut = path.join(outdir, 'theme');
    const themeOptions = {
      ...options,
      entryPoints: [themeEntry],
      outfile: path.join(themeOut, `${theme}.js`),
      globalName: undefined,
      plugins: [
        sourceJsExtensionPlugin,
        sassPlugin({ type: 'css', embedded: false, importers: [sassImporter], silenceDeprecations: ['import'] }),
      ],
      sourcemap: false,
    };
    if (dev) {
      const ctx = await esbuild.context(themeOptions);
      await ctx.watch();
      contexts.push(ctx);
    } else await esbuild.build(themeOptions);
  }
  console.log(`${dev ? 'Watching' : 'Built'} ${pkg.name} → ${path.relative(root, outdir)}`);
}
const contexts = [];
await Promise.all(selected.map(buildPackage));
if (serve) {
  const mime = {
    '.html': 'text/html',
    '.js': 'text/javascript',
    '.css': 'text/css',
    '.json': 'application/json',
    '.png': 'image/png',
    '.jpg': 'image/jpeg',
    '.svg': 'image/svg+xml',
    '.ico': 'image/x-icon',
  };
  const certfile = path.join(root, 'dev.pem');
  const keyfile = path.join(root, 'dev-key.pem');
  const sslEnabled = fs.existsSync(certfile) && fs.existsSync(keyfile);
  const protocol = sslEnabled ? 'https' : 'http';
  const handler = async (req, res) => {
    try {
      const url = new URL(req.url, `${protocol}://localhost`);
      const pathname = decodeURIComponent(url.pathname);
      const resolved = path.resolve(root, `.${pathname}`);
      if (!resolved.startsWith(root + path.sep) && resolved !== root) {
        res.writeHead(403).end();
        return;
      }
      const stat = await fs.promises.stat(resolved);
      if (stat.isDirectory() && !url.pathname.endsWith('/')) {
        res.writeHead(302, { location: `${url.pathname}/${url.search}` }).end();
        return;
      }
      const file = stat.isDirectory() ? path.join(resolved, 'index.html') : resolved;
      const body = await fs.promises.readFile(file);
      res.writeHead(200, { 'content-type': mime[path.extname(file)] ?? 'application/octet-stream' }).end(body);
    } catch {
      res.writeHead(404).end('Not found');
    }
  };
  const server = sslEnabled
    ? https.createServer({ cert: fs.readFileSync(certfile), key: fs.readFileSync(keyfile) }, handler)
    : http.createServer(handler);
  await new Promise((resolve) => server.listen(port, '0.0.0.0', resolve));
  for (const pkg of selected) console.log(`${pkg.name}: ${protocol}://localhost:${port}/${pkg.path}/build/`);
  await new Promise((resolve) => {
    const close = async () => {
      server.close();
      await Promise.all(contexts.map((c) => c.dispose()));
      resolve();
    };
    process.once('SIGINT', close);
    process.once('SIGTERM', close);
  });
}
