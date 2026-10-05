/** Refresh the compact export/dependency inventory after ownership changes. */
import fs from 'node:fs';
import path from 'node:path';
import ts from 'typescript';
import { compactExports } from './workspace/exports.mjs';

const inventory = JSON.parse(fs.readFileSync('scripts/workspace/inventory.json', 'utf8'));
const packages = inventory.packages;
const byName = new Map(packages.map(p => [p.name, p]));
const catalog = new Set([...fs.readFileSync('pnpm-workspace.yaml', 'utf8').matchAll(/^  '([^']+)':/gm)].map(m => m[1]));
function walk(dir) {
    return fs.readdirSync(dir, { withFileTypes: true }).flatMap(e => e.isDirectory() ? walk(path.join(dir, e.name)) : [path.join(dir, e.name)]);
}
for (const pkg of packages) {
    const manifestPath = path.join(pkg.path, 'package.json');
    const manifest = JSON.parse(fs.readFileSync(manifestPath, 'utf8'));
    if (pkg.private !== true) {
        if (Array.isArray(manifest.files) && !manifest.files.includes('!**/*.tsbuildinfo')) {
            manifest.files.push('!**/*.tsbuildinfo');
            fs.writeFileSync(manifestPath, JSON.stringify(manifest, null, 2) + '\n');
        }
    }
    if (!fs.existsSync(`${pkg.path}/src`) || pkg.kind === 'distribution') continue;
    const file = manifestPath;
    const sources = walk(`${pkg.path}/src`).filter(s => !s.includes('/_test/'));
    const exports = {}, imports = new Set();
    for (const source of sources) {
        const leaf = source.slice((pkg.path + '/src/').length);
        if (/\.tsx?$/.test(source)) {
            const text = fs.readFileSync(source, 'utf8');
            const tree = ts.createSourceFile(source, text, ts.ScriptTarget.Latest, true);
            function visit(n) {
                if (ts.isStringLiteral(n) && ((ts.isImportDeclaration(n.parent) || ts.isExportDeclaration(n.parent)) && n.parent.moduleSpecifier === n || ts.isCallExpression(n.parent) && ['require', 'import'].includes(n.parent.expression.getText(tree)))) imports.add(n.text);
                if (ts.isImportTypeNode(n) && ts.isLiteralTypeNode(n.argument) && ts.isStringLiteral(n.argument.literal)) imports.add(n.argument.literal.text);
                ts.forEachChild(n, visit);
            }
            visit(tree);
            const subpath = leaf.replace(/\.(d\.ts|tsx?)$/, '');
            const declaration = leaf.endsWith('.d.ts');
            const entry = declaration ? { types: './src/' + leaf } : { 'molstar-src': './src/' + leaf, types: './lib/' + subpath + '.d.ts', import: './lib/' + subpath + '.js' };
            exports['./' + subpath] = entry;
            if (subpath.endsWith('/index')) exports['./' + subpath.slice(0, -6)] = entry;
            if (subpath === 'index') exports['.'] = entry;
        } else if (leaf.endsWith('.scss')) {
            exports['./' + leaf] = { 'molstar-src': './src/' + leaf, sass: './src/' + leaf, default: './lib/' + leaf };
            if (!path.basename(leaf).startsWith('_')) exports['./' + leaf.replace(/\.scss$/, '.css')] = './lib/' + leaf.replace(/\.scss$/, '.css');
        }
    }
    if (pkg.name === '@molstar/core') exports['./task'] = exports['./task/index'];
    if (pkg.name === '@molstar/plugin') exports['.'] = exports['./context'];
    if (pkg.name === '@molstar/viewer') exports['.'] = exports['./app'];
    if (pkg.name === '@molstar/mvs-builder') exports['.'] = exports['./mvs-data'];
    if (pkg.name === '@molstar/plugin-headless') exports['.'] = exports['./context'];
    const deps = Object.fromEntries(Object.entries(manifest.dependencies ?? {}).filter(([n]) => !byName.has(n)));
    for (const spec of imports) {
        if (spec.startsWith('.') || spec.startsWith('node:')) continue;
        const name = spec.startsWith('@') ? spec.split('/').slice(0, 2).join('/') : spec.split('/')[0];
        if (name === pkg.name) continue;
        if (byName.has(name)) deps[name] = 'workspace:*';
        else if (catalog.has(name) && !manifest.peerDependencies?.[name]) deps[name] = 'catalog:';
    }
    if (pkg.name === '@molstar/common-server') {
        deps['@types/serve-static'] = 'catalog:';
        deps['@types/express-serve-static-core'] = 'catalog:';
    }
    manifest.dependencies = deps;
    manifest.exports = compactExports(exports, walk(`${pkg.path}/src`).map(file => path.relative(`${pkg.path}/src`, file).split(path.sep).join('/')));
    fs.writeFileSync(file, JSON.stringify(manifest, null, 2) + '\n');
    const configPath = `${pkg.path}/tsconfig.json`;
    const config = JSON.parse(fs.readFileSync(configPath, 'utf8'));
    config.references = Object.keys(deps).filter(n => byName.has(n)).map(n => ({ path: path.relative(pkg.path, byName.get(n).path) }));
    fs.writeFileSync(configPath, JSON.stringify(config, null, 2) + '\n');
}
const solution = { files: [], references: packages.filter(p => fs.existsSync(`${p.path}/tsconfig.json`)).map(p => ({ path: './' + p.path })) };
fs.writeFileSync('tsconfig.json', JSON.stringify(solution, null, 2) + '\n');
console.log('Refreshed workspace exports and project references.');
