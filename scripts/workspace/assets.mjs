import fs from 'node:fs/promises';
import path from 'node:path';
import { fileURLToPath } from 'node:url';
import * as sass from 'sass';

const root = path.resolve(path.dirname(fileURLToPath(import.meta.url)), '../..');
const inventory = JSON.parse(await fs.readFile(path.join(root, 'scripts/workspace/inventory.json'), 'utf8'));
const sassImporter = new sass.NodePackageImporter(root);
const selected = process.argv.slice(2).filter(x => !x.startsWith('--'));
const packages = (inventory.packages ?? []).filter(p => !selected.length || selected.includes(p.name) || selected.includes(p.path));
if (selected.length && packages.length !== selected.length) throw new Error(`Unknown package(s): ${selected.filter(name => !packages.some(pkg => pkg.name === name || pkg.path === name)).join(', ')}`);
const ignored = /(^|\/)(node_modules|lib|dist|build|_test|__tests__|tests?|fixtures?|\.git)(\/|$)|\.test\.[^.]+$/i;
const assetExtensions = new Set(['.scss', '.sass', '.css', '.html', '.htm', '.json', '.wasm', '.glsl', '.vert', '.frag', '.png', '.jpg', '.jpeg', '.gif', '.svg', '.ico', '.woff', '.woff2', '.ttf', '.otf', '.bin', '.dat', '.csv', '.txt', '.xml']);
async function walk(dir, base, found = []) {
    for (const entry of await fs.readdir(dir, { withFileTypes: true })) {
        const abs = path.join(dir, entry.name), rel = path.relative(base, abs).split(path.sep).join('/');
        if (ignored.test(rel)) continue;
        if (entry.isDirectory()) await walk(abs, base, found);
        else if (assetExtensions.has(path.extname(entry.name).toLowerCase())) found.push({ abs, rel });
    }
    return found;
}
for (const pkg of packages) {
    const source = path.join(root, pkg.path, 'src');
    const output = path.join(root, pkg.path, 'lib');
    try { await fs.access(source); } catch { continue; }
    for (const { abs, rel } of await walk(source, source)) {
        const isSass = /\.s[ac]ss$/i.test(abs);
        const target = path.join(output, rel);
        await fs.mkdir(path.dirname(target), { recursive: true });
        if (isSass && !path.basename(abs).startsWith('_')) {
            await fs.copyFile(abs, target);
            const result = await sass.compileAsync(abs, { style: 'expanded', loadPaths: [source], importers: [sassImporter] });
            await fs.writeFile(target.replace(/\.s[ac]ss$/i, '.css'), result.css);
        } else await fs.copyFile(abs, target);
    }
}
