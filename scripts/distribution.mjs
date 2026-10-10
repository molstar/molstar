import * as esbuild from 'esbuild';
import fs from 'node:fs/promises';
import fsSync from 'node:fs';
import path from 'node:path';
import { fileURLToPath, pathToFileURL } from 'node:url';
import * as sass from 'sass';
import { sassPlugin } from 'esbuild-sass-plugin';
import { resolveExport as resolveExportMap } from './workspace/exports.mjs';

const root = path.resolve(path.dirname(fileURLToPath(import.meta.url)), '..');
const inventory = JSON.parse(await fs.readFile(path.join(root, 'scripts/workspace/inventory.json'), 'utf8'));
const sassImporter = new sass.NodePackageImporter(root);
const packages = new Map((inventory.packages ?? []).map((p) => [p.name, p]));
const sassPluginImporter = {
  canonicalize(specifier, { containingUrl }) {
    let file;
    if (specifier.startsWith('pkg:')) {
      const request = specifier.slice(4);
      const separator = request.startsWith('@') ? request.indexOf('/', request.indexOf('/') + 1) : request.indexOf('/');
      if (separator < 0) return null;
      const pkg = packages.get(request.slice(0, separator));
      const subpath = request.slice(separator + 1);
      if (!pkg || !subpath) return null;
      const sourceRoot = path.resolve(root, pkg.path, 'src');
      file = path.resolve(sourceRoot, subpath);
      if (!file.startsWith(sourceRoot + path.sep)) return null;
    } else if (containingUrl?.protocol === 'file:') {
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
    const resolved = candidates.find((candidate) => fsSync.existsSync(candidate));
    return resolved ? pathToFileURL(resolved) : null;
  },
  load(url) {
    const file = fileURLToPath(url);
    let contents = fsSync.readFileSync(file, 'utf8');
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
        const resolved = candidates.find((candidate) => fsSync.existsSync(candidate));
        return resolved ? `${prefix}${quote}${pathToFileURL(resolved).href}${quote}` : match;
      },
    );
    return { contents, syntax: file.endsWith('.sass') ? 'indented' : 'scss' };
  },
};
const destination = path.join(root, 'distributions/molstar/build');
const esmDestination = path.join(destination, 'esm');
const baseArg = process.argv.indexOf('--base-url');
const baseUrlArg = baseArg >= 0 ? process.argv[baseArg + 1] : '/build/esm/';
if (!baseUrlArg || (!baseUrlArg.startsWith('/') && !/^https?:\/\//u.test(baseUrlArg)))
  throw new Error('--base-url must be an absolute path or HTTP(S) URL');
const baseUrl = `${baseUrlArg.replace(/\/+$/u, '')}/`;
const appDirs = [
  { package: 'apps/viewer', target: 'viewer' },
  { package: 'apps/mvs-stories', target: 'mvs-stories' },
];
// Public browser entrypoints kept intentionally small; aliases are exact import-map keys.
const supported = [
  {
    alias: '@molstar/model/model/structure',
    pkg: '@molstar/model',
    subpath: './model/structure',
    output: 'model/model/structure',
  },
  { alias: '@molstar/core/state', pkg: '@molstar/core', subpath: './state', output: 'core/state' },
  { alias: '@molstar/plugin/config', pkg: '@molstar/plugin', subpath: './config', output: 'plugin/config' },
  {
    alias: '@molstar/plugin/state/formats/trajectory/pdb',
    pkg: '@molstar/plugin',
    subpath: './state/formats/trajectory/pdb',
    output: 'plugin/state/formats/trajectory/pdb',
  },
  {
    alias: '@molstar/plugin/state/transforms/catalog',
    pkg: '@molstar/plugin',
    subpath: './state/transforms/catalog',
    output: 'plugin/state/transforms/catalog',
  },
  { alias: '@molstar/plugin-ui/context', pkg: '@molstar/plugin-ui', subpath: './context', output: 'plugin-ui/context' },
  { alias: '@molstar/core/util/color', pkg: '@molstar/core', subpath: './util/color', output: 'core/util/color' },
  { alias: '@molstar/core/task', pkg: '@molstar/core', subpath: './task', output: 'core/task' },
  { alias: '@molstar/io/reader/cif', pkg: '@molstar/io', subpath: './reader/cif', output: 'io/reader/cif' },
  { alias: '@molstar/plugin', pkg: '@molstar/plugin', subpath: '.', output: 'plugin/index' },
  { alias: '@molstar/plugin/spec', pkg: '@molstar/plugin', subpath: './spec', output: 'plugin/spec' },
  {
    alias: '@molstar/plugin/default-spec',
    pkg: '@molstar/plugin',
    subpath: './default-spec',
    output: 'plugin/default-spec',
  },
  { alias: '@molstar/plugin-ui', pkg: '@molstar/plugin-ui', subpath: '.', output: 'plugin-ui/index' },
  { alias: '@molstar/plugin-ui/react18', pkg: '@molstar/plugin-ui', subpath: './react18', output: 'plugin-ui/react18' },
  { alias: '@molstar/plugin-ui/spec', pkg: '@molstar/plugin-ui', subpath: './spec', output: 'plugin-ui/spec' },
  {
    alias: '@molstar/plugin-ui/default-spec',
    pkg: '@molstar/plugin-ui',
    subpath: './default-spec',
    output: 'plugin-ui/default-spec',
  },
];

async function exists(p) {
  try {
    await fs.access(p);
    return true;
  } catch {
    return false;
  }
}
async function copyTree(source, target) {
  if (!(await exists(source)))
    throw new Error(`Missing app build output ${path.relative(root, source)}. Build app packages first.`);
  await fs.mkdir(target, { recursive: true });
  for (const entry of await fs.readdir(source, { withFileTypes: true })) {
    const from = path.join(source, entry.name),
      to = path.join(target, entry.name);
    if (entry.isDirectory()) await copyTree(from, to);
    else await fs.copyFile(from, to);
  }
}
function getSourceTarget(value) {
  if (typeof value === 'string') return undefined;
  if (!value || typeof value !== 'object') return undefined;
  if (typeof value['molstar-src'] === 'string') return value['molstar-src'];
  for (const v of Object.values(value)) {
    const result = getSourceTarget(v);
    if (result) return result;
  }
}
function resolveExport(pkgName, subpath) {
  const entry = packages.get(pkgName);
  if (!entry) throw new Error(`Inventory has no package ${pkgName}`);
  const dir = path.join(root, entry.path);
  const manifest = JSON.parse(fsSync.readFileSync(path.join(dir, 'package.json'), 'utf8'));
  const exports = manifest.exports ?? {};
  const key = subpath === '.' ? '.' : subpath;
  const record = resolveExportMap(exports, key);
  if (!record) throw new Error(`${pkgName} does not export ${key}`);
  const target = getSourceTarget(record);
  if (!target) throw new Error(`${pkgName} ${key} has no molstar-src export condition`);
  return path.resolve(dir, target);
}
const sourceJsExtensionPlugin = {
  name: 'molstar-source-js-extension',
  setup(build) {
    build.onResolve({ filter: /^\./ }, (args) => {
      if (!args.path.endsWith('.js')) return null;
      const file = path.resolve(args.resolveDir, args.path);
      if (fsSync.existsSync(file)) return null;
      for (const ext of ['.ts', '.tsx']) {
        const source = file.slice(0, -'.js'.length) + ext;
        if (fsSync.existsSync(source)) return { path: source };
      }
      return null;
    });
  },
};
// Clear only generated distribution output.
await fs.rm(destination, { recursive: true, force: true });
await fs.mkdir(destination, { recursive: true });
for (const app of appDirs) {
  const source = path.join(root, app.package, 'build');
  await copyTree(source, path.join(destination, app.target));
}

const viewerEntry = path.join(root, 'apps/viewer/src/app.ts');
const viewerEntryTsx = viewerEntry + 'x';
const viewerSource = (await exists(viewerEntry)) ? viewerEntry : viewerEntryTsx;
if (!(await exists(viewerSource))) throw new Error('Missing apps/viewer/src/app.ts(x) modular ESM entry');
const entryPoints = [{ in: viewerSource, out: 'viewer' }];
const importMap = { '@molstar/viewer': `${baseUrl}viewer.js` };
for (const entry of supported) {
  const source = resolveExport(entry.pkg, entry.subpath);
  entryPoints.push({ in: source, out: entry.output });
  importMap[entry.alias] = `${baseUrl}${entry.output}.js`;
}
await esbuild.build({
  absWorkingDir: root,
  entryPoints,
  outdir: esmDestination,
  bundle: true,
  splitting: true,
  format: 'esm',
  platform: 'browser',
  target: ['es2022'],
  conditions: ['molstar-src', 'import', 'default'],
  chunkNames: 'chunks/[name]-[hash]',
  assetNames: 'assets/[name]-[hash]',
  tsconfig: path.join(root, 'tsconfig.base.json'),
  minify: true,
  minifyIdentifiers: false,
  sourcemap: false,
  external: ['crypto', 'fs', 'path', 'stream'],
  loader: {
    '.png': 'file',
    '.jpg': 'file',
    '.jpeg': 'file',
    '.gif': 'file',
    '.svg': 'file',
    '.woff': 'file',
    '.woff2': 'file',
    '.ttf': 'file',
    '.otf': 'file',
    '.ico': 'file',
  },
  plugins: [
    sourceJsExtensionPlugin,
    sassPlugin({ type: 'css', embedded: false, importers: [sassPluginImporter], silenceDeprecations: ['import'] }),
  ],
  define: {
    'process.env.NODE_ENV': '"production"',
    __MOLSTAR_PLUGIN_VERSION__: JSON.stringify(
      JSON.parse(await fs.readFile(path.join(root, 'version.json'), 'utf8')).version,
    ),
  },
  logLevel: 'info',
});
await fs.writeFile(
  path.join(esmDestination, 'import-map.json'),
  JSON.stringify({ imports: importMap }, null, 2) + '\n',
);

// Emit the primary plugin skin as a distributable stylesheet when it is available.
const ui = packages.get('@molstar/plugin-ui');
if (ui) {
  const uiRoot = path.join(root, ui.path);
  const candidates = [path.join(uiRoot, 'src/skin/light.scss'), path.join(uiRoot, 'src/skin/light.sass')];
  let skinPath;
  for (const candidate of candidates)
    if (await exists(candidate)) {
      skinPath = candidate;
      break;
    }
  if (skinPath) {
    const result = await sass.compileAsync(skinPath, {
      style: 'compressed',
      loadPaths: [path.join(uiRoot, 'src')],
      importers: [sassImporter],
    });
    await fs.writeFile(path.join(esmDestination, 'viewer.css'), result.css);
  }
}
const classicViewerCss = path.join(root, 'apps/viewer/build/molstar.css');
if (!(await exists(path.join(esmDestination, 'viewer.css'))) && (await exists(classicViewerCss))) {
  await fs.copyFile(classicViewerCss, path.join(esmDestination, 'viewer.css'));
}
console.log(`Staged Viewer and MVS Stories plus browser ESM into ${path.relative(root, destination)}.`);
