import * as esbuild from 'esbuild';
import fs from 'node:fs';
import path from 'node:path';
import http from 'node:http';
import { fileURLToPath, pathToFileURL } from 'node:url';
import { sassPlugin } from 'esbuild-sass-plugin';

const root = path.resolve(path.dirname(fileURLToPath(import.meta.url)), '../..');
const args = process.argv.slice(2);
const value = (flag, fallback) => { const i = args.indexOf(flag); return i >= 0 ? args[i + 1] : fallback; };
const all = args.includes('--all');
const dev = args.includes('--dev') || args.includes('--watch');
const production = args.includes('--prd') || !dev;
const serve = dev || args.includes('--serve');
const port = Number(value('--port', '1338'));
const target = args.find(a => !a.startsWith('-'));
const timestamp = Number(process.env.MOLSTAR_BUILD_TIMESTAMP ?? Date.now());
const version = JSON.parse(fs.readFileSync(path.join(root, 'version.json'), 'utf8')).version;
const inventoryPath = path.join(root, 'scripts/workspace/inventory.json');
const inventory = fs.existsSync(inventoryPath) ? JSON.parse(fs.readFileSync(inventoryPath, 'utf8')) : { packages: [] };
const packagesByName = new Map((inventory.packages ?? []).map(p => [p.name, p]));
const sassImporter = {
    canonicalize(specifier, { containingUrl }) {
        let file;
        if (specifier.startsWith('pkg:')) {
            const request = specifier.slice(4);
            const packageName = request.startsWith('@') ? request.indexOf('/', request.indexOf('/') + 1) : request.indexOf('/');
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
        const candidates = ext ? [file] : [file + '.scss', path.join(path.dirname(file), `_${path.basename(file)}.scss`), file + '.sass', path.join(path.dirname(file), `_${path.basename(file)}.sass`)];
        const resolved = candidates.find(candidate => fs.existsSync(candidate));
        return resolved ? pathToFileURL(resolved) : null;
    },
    load(url) {
        const file = fileURLToPath(url);
        let contents = fs.readFileSync(file, 'utf8');
        contents = contents.replace(/(@(?:use|forward|import)\s+)(['"])(\.{1,2}\/[^'"]+)\2/gu, (match, prefix, quote, specifier) => {
            const base = path.resolve(path.dirname(file), specifier);
            const ext = path.extname(base);
            const candidates = ext ? [base, path.join(path.dirname(base), `_${path.basename(base)}`)] : [base + '.scss', path.join(path.dirname(base), `_${path.basename(base)}.scss`), base + '.sass', path.join(path.dirname(base), `_${path.basename(base)}.sass`)];
            const resolved = candidates.find(candidate => fs.existsSync(candidate));
            return resolved ? `${prefix}${quote}${pathToFileURL(resolved).href}${quote}` : match;
        });
        return { contents, syntax: file.endsWith('.sass') ? 'indented' : 'scss' };
    }
};
const pkgInfo = p => {
    const manifest = JSON.parse(fs.readFileSync(path.join(root, p.path, 'package.json'), 'utf8'));
    return { ...p, ...manifest, config: manifest.molstar ?? {} };
};
const apps = (inventory.packages ?? []).filter(p => {
    if (!['app', 'example'].includes(p.kind)) return false;
    const manifestPath = path.join(root, p.path, 'package.json');
    if (!fs.existsSync(manifestPath)) return false;
    const manifest = JSON.parse(fs.readFileSync(manifestPath, 'utf8'));
    return manifest.molstar?.platform !== 'node';
});
const selected = all ? apps : target ? apps.filter(p => p.name === target || p.path === target) : [];
if (!selected.length) {
    console.error('Usage: node scripts/esbuild/app.mjs <package-name|path> [--prd|--dev] | --all [--prd|--dev]');
    process.exit(2);
}
for (const p of selected) if (!fs.existsSync(path.join(root, p.path, 'package.json'))) throw new Error(`Missing manifest: ${p.path}/package.json`);

function staticAssetsPlugin(outdir) {
    return {
        name: 'molstar-static-assets',
        setup(build) {
            build.onLoad({ filter: /\.(jpg|jpeg|png|gif|ico|html|htm|svg|wasm|bin|dat)$/i }, async ({ path: input }) => {
                const ext = path.extname(input).toLowerCase();
                const name = path.basename(input);
                const isImage = /\.(jpg|jpeg|png|gif|svg)$/i.test(ext);
                const destDir = isImage ? path.join(outdir, 'images') : outdir;
                await fs.promises.mkdir(destDir, { recursive: true });
                await fs.promises.copyFile(input, path.join(destDir, name));
                if (/\.(html|htm|ico)$/i.test(ext)) return { contents: '', loader: 'empty' };
                return { contents: `${isImage ? 'images/' : ''}${name}`, loader: 'text' };
            });
        }
    };
}
const sourceJsExtensionPlugin = {
    name: 'molstar-source-js-extension',
    setup(build) {
        build.onResolve({ filter: /^\./ }, args => {
            if (!args.path.endsWith('.js')) return null;
            const file = path.resolve(args.resolveDir, args.path);
            if (fs.existsSync(file)) return null;
            for (const ext of ['.ts', '.tsx']) {
                const source = file.slice(0, -'.js'.length) + ext;
                if (fs.existsSync(source)) return { path: source };
            }
            return null;
        });
    }
};
function entryFor(pkg, config) {
    const raw = path.resolve(root, pkg.path, config.entry ?? 'src/index.ts');
    if (fs.existsSync(raw)) return raw;
    if (fs.existsSync(raw + 'x')) return raw + 'x';
    throw new Error(`${pkg.name}: configured entry does not exist: ${path.relative(root, raw)}`);
}
async function buildPackage(rawPkg) {
    const pkg = pkgInfo(rawPkg), cfg = pkg.config;
    const dir = path.join(root, pkg.path), outdir = path.join(dir, 'build');
    const entry = entryFor(pkg, cfg);
    await fs.promises.mkdir(outdir, { recursive: true });
    const filename = cfg.filename ?? (pkg.kind === 'example' ? 'index.js' : 'molstar.js');
    const outfile = path.join(outdir, filename);
    const options = {
        absWorkingDir: root,
        entryPoints: [entry], bundle: true, platform: 'browser', format: 'iife',
        globalName: cfg.globalName ?? 'molstar', outfile,
        tsconfig: fs.existsSync(path.join(dir, 'tsconfig.json')) ? path.join(dir, 'tsconfig.json') : path.join(root, 'tsconfig.json'),
        minify: production, minifyIdentifiers: false, sourcemap: !production,
        conditions: ['molstar-src', 'import', 'default'],
        external: ['crypto', 'fs', 'path', 'stream'],
        plugins: [sourceJsExtensionPlugin, staticAssetsPlugin(outdir), sassPlugin({ type: 'css', embedded: false, importers: [sassImporter], silenceDeprecations: ['import'] })],
        define: {
            'process.env.NODE_ENV': JSON.stringify(production ? 'production' : 'development'),
            'process.env.DEBUG': JSON.stringify(process.env.DEBUG || false),
            __MOLSTAR_PLUGIN_VERSION__: JSON.stringify(version),
            __MOLSTAR_BUILD_TIMESTAMP__: String(timestamp),
        }, logLevel: 'info'
    };
    if (dev) {
        const ctx = await esbuild.context(options);
        await ctx.watch();
        contexts.push(ctx);
    } else await esbuild.build(options);
    for (const theme of cfg.themes ?? []) {
        const themeRoot = path.join(dir, 'src/theme');
        let themeEntry = path.join(themeRoot, `${theme}.ts`);
        if (!fs.existsSync(themeEntry)) themeEntry += 'x';
        if (!fs.existsSync(themeEntry)) throw new Error(`${pkg.name}: theme entry not found for ${theme}`);
        const themeOut = path.join(outdir, 'theme');
        const themeOptions = { ...options, entryPoints: [themeEntry], outfile: path.join(themeOut, `${theme}.js`), globalName: undefined, plugins: [sourceJsExtensionPlugin, sassPlugin({ type: 'css', embedded: false, importers: [sassImporter], silenceDeprecations: ['import'] })], sourcemap: false };
        if (dev) { const ctx = await esbuild.context(themeOptions); await ctx.watch(); contexts.push(ctx); } else await esbuild.build(themeOptions);
    }
    if (pkg.kind === 'example') {
        const indexCss = path.join(outdir, 'index.css');
        if (fs.existsSync(indexCss)) await fs.promises.rename(indexCss, path.join(outdir, 'molstar.css'));
    }
    console.log(`${dev ? 'Watching' : 'Built'} ${pkg.name} → ${path.relative(root, outdir)}`);
}
const contexts = [];
await Promise.all(selected.map(buildPackage));
if (serve) {
    const mime = { '.html': 'text/html', '.js': 'text/javascript', '.css': 'text/css', '.json': 'application/json', '.png': 'image/png', '.jpg': 'image/jpeg', '.svg': 'image/svg+xml', '.ico': 'image/x-icon' };
    const server = http.createServer(async (req, res) => {
        try {
            const pathname = decodeURIComponent(new URL(req.url, 'http://localhost').pathname);
            const resolved = path.resolve(root, `.${pathname}`);
            if (!resolved.startsWith(root + path.sep) && resolved !== root) { res.writeHead(403).end(); return; }
            const stat = await fs.promises.stat(resolved);
            const file = stat.isDirectory() ? path.join(resolved, 'index.html') : resolved;
            const body = await fs.promises.readFile(file);
            res.writeHead(200, { 'content-type': mime[path.extname(file)] ?? 'application/octet-stream' }).end(body);
        } catch { res.writeHead(404).end('Not found'); }
    });
    await new Promise(resolve => server.listen(port, '0.0.0.0', resolve));
    console.log(`Static server: http://localhost:${port}`);
    await new Promise(resolve => {
        const close = async () => { server.close(); await Promise.all(contexts.map(c => c.dispose())); resolve(); };
        process.once('SIGINT', close); process.once('SIGTERM', close);
    });
}
