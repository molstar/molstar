import fs from 'node:fs';
import path from 'node:path';
import { fileURLToPath } from 'node:url';
import { builtinModules } from 'node:module';
import ts from 'typescript';
import { expandExports, exportTargets } from './exports.mjs';

const root = path.resolve(path.dirname(fileURLToPath(import.meta.url)), '../..');
const inventory = JSON.parse(fs.readFileSync(path.join(root, 'scripts/workspace/inventory.json'), 'utf8'));
const version = JSON.parse(fs.readFileSync(path.join(root, 'version.json'), 'utf8')).version;
const sourceOnly = process.argv.includes('--source-only');
const packages = inventory.packages ?? [];
const publicPackages = packages.filter(p => p.private !== true);
const byName = new Map(packages.map(p => [p.name, p]));
const errors = [];
const fields = ['dependencies', 'devDependencies', 'peerDependencies', 'optionalDependencies'];
const manifestCache = new Map();
if (!Array.isArray(inventory.packages)) errors.push('inventory.packages must be an array');
if (byName.size !== packages.length) errors.push('inventory contains duplicate package names');
if (new Set(packages.map(pkg => pkg.path)).size !== packages.length) errors.push('inventory contains duplicate package paths');
for (const pkg of packages) {
    if (!pkg.name || !pkg.path || path.isAbsolute(pkg.path) || pkg.path.split(/[\\/]/).includes('..')) errors.push(`invalid package inventory entry: ${JSON.stringify(pkg)}`);
    const manifest = manifestFor(pkg);
    if (manifest && manifest.name !== pkg.name) errors.push(`${pkg.path}: manifest name ${manifest.name} != inventory name ${pkg.name}`);
}
function packageManifestPaths(dir, results = []) {
    for (const ent of fs.readdirSync(dir, { withFileTypes: true })) {
        if (ent.isDirectory() && ['node_modules', '.git', '.pnpm-store', '.cache', 'build', 'lib', '.v6'].includes(ent.name)) continue;
        const abs = path.join(dir, ent.name);
        if (ent.isDirectory()) packageManifestPaths(abs, results);
        else if (ent.isFile() && ent.name === 'package.json' && abs !== path.join(root, 'package.json')) results.push(path.relative(root, path.dirname(abs)).split(path.sep).join('/'));
    }
    return results;
}
const inventoryPaths = new Set(packages.map(pkg => pkg.path));
for (const manifestDir of packageManifestPaths(root)) if (!inventoryPaths.has(manifestDir)) errors.push(`package manifest missing from inventory: ${manifestDir}/package.json`);
function manifestFor(pkg) {
    if (!manifestCache.has(pkg.name)) {
        const file = path.join(root, pkg.path, 'package.json');
        manifestCache.set(pkg.name, fs.existsSync(file) ? JSON.parse(fs.readFileSync(file, 'utf8')) : undefined);
    }
    return manifestCache.get(pkg.name);
}
function resolveFile(file) {
    if (fs.existsSync(file) && fs.statSync(file).isFile()) return true;
    for (const ext of ['.js', '.mjs', '.cjs', '.d.ts', '.ts', '.tsx', '.json', '.css', '.scss']) {
        if (fs.existsSync(file + ext)) return true;
    }
    return ['index.js', 'index.d.ts', 'index.ts', 'index.tsx', 'index.css'].some(name => fs.existsSync(path.join(file, name)));
}
function sourceCounterpart(pkg, target) {
    if (!target.startsWith('./lib/')) return undefined;
    const sourceBase = path.join(root, pkg.path, target.replace(/^\.\/lib\//u, 'src/'));
    const candidates = [sourceBase];
    for (const [from, to] of [['.js', '.ts'], ['.js', '.tsx'], ['.mjs', '.mts'], ['.cjs', '.cts'], ['.d.ts', '.ts'], ['.d.ts', '.tsx'], ['.css', '.scss'], ['.css', '.sass']]) {
        if (sourceBase.endsWith(from)) candidates.push(sourceBase.slice(0, -from.length) + to);
    }
    return candidates.find(resolveFile);
}
function listFiles(dir, collected = []) {
    if (!fs.existsSync(dir)) return collected;
    for (const ent of fs.readdirSync(dir, { withFileTypes: true })) {
        if (['node_modules', '.git'].includes(ent.name)) continue;
        const abs = path.join(dir, ent.name);
        if (ent.isDirectory()) listFiles(abs, collected);
        else collected.push(abs);
    }
    return collected;
}
function validateExport(pkg, key, target) {
    if (!target.startsWith('./') || target.split('/').includes('..')) {
        errors.push(`${pkg.name}: export ${key} target must stay package-local: ${target}`);
        return;
    }
    const dir = path.join(root, pkg.path);
    const output = path.resolve(dir, target);
    if (!resolveFile(output) && !(sourceOnly && sourceCounterpart(pkg, target))) errors.push(`${pkg.name}: export ${key} target does not resolve: ${target}`);
}
for (const pkg of publicPackages) {
    const manifest = manifestFor(pkg);
    if (!manifest) { errors.push(`${pkg.name}: missing manifest at ${pkg.path}/package.json`); continue; }
    if (manifest.name !== pkg.name) errors.push(`${pkg.path}: manifest name ${manifest.name} != ${pkg.name}`);
    if (manifest.version !== version) errors.push(`${pkg.name}: version ${manifest.version} != ${version}`);
    if (pkg.kind !== 'distribution') {
        if (!manifest.exports || typeof manifest.exports !== 'object' || !Object.keys(manifest.exports).length) errors.push(`${pkg.name}: public package has no exports map`);
        try {
            const dir = path.join(root, pkg.path);
            const files = listFiles(path.join(dir, 'src')).concat(sourceOnly ? [] : listFiles(path.join(dir, 'lib'))).map(file => path.relative(dir, file).split(path.sep).join('/'));
            if (sourceOnly) for (const file of [...files]) {
                if (/^src\/.*\.tsx?$/.test(file) && !file.endsWith('.d.ts')) {
                    const base = file.replace(/^src\//u, 'lib/').replace(/\.tsx?$/u, '');
                    files.push(`${base}.js`, `${base}.d.ts`);
                } else if (/^src\/.*\.s[ac]ss$/.test(file)) {
                    files.push(file.replace(/^src\//u, 'lib/'));
                    if (!path.basename(file).startsWith('_')) files.push(file.replace(/^src\//u, 'lib/').replace(/\.s[ac]ss$/u, '.css'));
                }
            }
            for (const [key, value] of Object.entries(expandExports(manifest.exports ?? {}, files))) {
                for (const target of exportTargets(value)) validateExport(pkg, key, target);
            }
        } catch (error) { errors.push(`${pkg.name}: ${error.message}`); }
    }
}
for (const pkg of packages) {
    const manifest = manifestFor(pkg);
    if (!manifest) continue;
    for (const field of fields) for (const [name, range] of Object.entries(manifest[field] ?? {})) {
        if (byName.has(name) && range !== 'workspace:*') errors.push(`${pkg.name}: ${field}.${name} must use workspace:* (found ${range})`);
    }
}

const sourceFiles = new Map();
function scan(dir, pkg) {
    if (!fs.existsSync(dir)) return;
    for (const ent of fs.readdirSync(dir, { withFileTypes: true })) {
        if (ent.isDirectory() && /^(node_modules|lib|dist|build|tests?|_test|__tests__|fixtures?)$/i.test(ent.name)) continue;
        const abs = path.join(dir, ent.name);
        if (ent.isDirectory()) scan(abs, pkg);
        else if (/\.[cm]?tsx?$/.test(ent.name) && !/\.test\.[cm]?tsx?$/.test(ent.name)) sourceFiles.set(abs, pkg);
    }
}
for (const pkg of packages) scan(path.join(root, pkg.path, 'src'), pkg);
const graph = new Map(packages.map(pkg => [pkg.name, new Set()]));
const sourceOwner = file => {
    let owner;
    for (const pkg of packages) {
        const dir = path.resolve(root, pkg.path);
        if (file === dir || file.startsWith(dir + path.sep)) {
            if (!owner || dir.length > path.resolve(root, owner.path).length) owner = pkg;
        }
    }
    return owner;
};
function packageName(specifier) {
    if (specifier.startsWith('@')) return specifier.split('/').slice(0, 2).join('/');
    return specifier.split('/')[0];
}
function ambientTypesPackage(name) {
    return name.startsWith('@') ? `@types/${name.slice(1).replace('/', '__')}` : `@types/${name}`;
}
function importsFrom(file, source) {
    const result = [];
    function isTypeOnlyImport(node) {
        const clause = node.importClause;
        if (clause?.isTypeOnly) return true;
        const names = [];
        if (clause?.name) names.push(clause.name.text);
        if (clause?.namedBindings && ts.isNamespaceImport(clause.namedBindings)) names.push(clause.namedBindings.name.text);
        if (clause?.namedBindings && ts.isNamedImports(clause.namedBindings)) names.push(...clause.namedBindings.elements.map(element => element.name.text));
        if (!names.length) return false;
        const uses = new Map(names.map(name => [name, []]));
        const findUses = current => {
            let declarationParent = current;
            while (declarationParent && declarationParent !== source && declarationParent !== node.importClause) declarationParent = declarationParent.parent;
            if (ts.isIdentifier(current) && uses.has(current.text) && declarationParent !== node.importClause) {
                let parent = current.parent, typePosition = false;
                while (parent && parent !== source) {
                    if (ts.isTypeNode(parent)) { typePosition = true; break; }
                    if (ts.isExpression(parent)) break;
                    parent = parent.parent;
                }
                uses.get(current.text).push(typePosition);
            }
            ts.forEachChild(current, findUses);
        };
        findUses(source);
        const found = [...uses.values()].flat();
        return found.length > 0 && found.every(Boolean);
    }
    const visit = node => {
        if (ts.isImportDeclaration(node) && node.moduleSpecifier && ts.isStringLiteral(node.moduleSpecifier)) result.push({ specifier: node.moduleSpecifier.text, typeOnly: isTypeOnlyImport(node) });
        if (ts.isExportDeclaration(node) && node.moduleSpecifier && ts.isStringLiteral(node.moduleSpecifier)) result.push({ specifier: node.moduleSpecifier.text, typeOnly: node.isTypeOnly });
        if (ts.isImportEqualsDeclaration(node) && ts.isExternalModuleReference(node.moduleReference) && node.moduleReference.expression && ts.isStringLiteral(node.moduleReference.expression)) result.push({ specifier: node.moduleReference.expression.text, typeOnly: false });
        if (ts.isCallExpression(node) && node.expression.kind === ts.SyntaxKind.ImportKeyword && node.arguments.length === 1 && ts.isStringLiteral(node.arguments[0])) result.push({ specifier: node.arguments[0].text, typeOnly: false });
        if (ts.isCallExpression(node) && ts.isIdentifier(node.expression) && node.expression.text === 'require' && node.arguments.length === 1 && ts.isStringLiteral(node.arguments[0])) result.push({ specifier: node.arguments[0].text, typeOnly: false });
        if (ts.isImportTypeNode(node) && ts.isLiteralTypeNode(node.argument) && ts.isStringLiteral(node.argument.literal)) result.push({ specifier: node.argument.literal.text, typeOnly: true });
        ts.forEachChild(node, visit);
    };
    visit(source);
    return result;
}
for (const [file, pkg] of sourceFiles) {
    const text = fs.readFileSync(file, 'utf8');
    const source = ts.createSourceFile(file, text, ts.ScriptTarget.Latest, true, file.endsWith('x') ? ts.ScriptKind.TSX : ts.ScriptKind.TS);
    const rel = path.relative(path.join(root, pkg.path), file).split(path.sep).join('/');
    for (const imported of importsFrom(file, source)) {
        const { specifier } = imported;
        if (specifier.startsWith('.')) {
            const target = path.resolve(path.dirname(file), specifier);
            const owner = sourceOwner(target);
            if (owner && owner.name !== pkg.name) errors.push(`${pkg.name}: cross-package relative import ${specifier} in ${rel} (resolves into ${owner.name})`);
            continue;
        }
        if (specifier.startsWith('#')) continue;
        const dependency = byName.get(packageName(specifier));
        if (dependency && dependency.name !== pkg.name) {
            const manifest = manifestFor(pkg);
            if (!manifest || !fields.some(field => Object.hasOwn(manifest[field] ?? {}, dependency.name))) errors.push(`${pkg.name}: missing direct dependency on ${dependency.name} (${specifier} in ${rel})`);
            graph.get(pkg.name)?.add(dependency.name);
        } else if (!dependency && !specifier.startsWith('node:') && !ts.isExternalModuleNameRelative(specifier)) {
            const name = packageName(specifier);
            if (builtinModules.includes(name) || builtinModules.includes(`node:${name}`) || name.startsWith('#')) continue;
            const manifest = manifestFor(pkg);
            if (name.startsWith('@molstar/')) errors.push(`${pkg.name}: imported workspace package ${name} is missing from inventory (${specifier} in ${rel})`);
            else if (!manifest || (!fields.some(field => Object.hasOwn(manifest[field] ?? {}, name)) && !(imported.typeOnly && fields.some(field => Object.hasOwn(manifest[field] ?? {}, ambientTypesPackage(name)))))) errors.push(`${pkg.name}: undeclared external import ${specifier} in ${rel}`);
        }
    }
}
const state = new Map(), stack = [];
function visit(name) {
    if (state.get(name) === 1) {
        errors.push(`package dependency cycle: ${stack.slice(stack.indexOf(name)).concat(name).join(' -> ')}`);
        return;
    }
    if (state.get(name) === 2) return;
    state.set(name, 1); stack.push(name);
    for (const child of graph.get(name) ?? []) visit(child);
    stack.pop(); state.set(name, 2);
}
for (const pkg of packages) visit(pkg.name);
if (errors.length) {
    const priority = error => /package dependency cycle|cross-package relative import|missing direct dependency|undeclared external import|missing from inventory/.test(error) ? 0 : 1;
    const ordered = [...errors].sort((a, b) => priority(a) - priority(b));
    const shown = ordered.slice(0, 150);
    console.error(shown.join('\n'));
    if (ordered.length > shown.length) console.error(`... ${ordered.length - shown.length} additional validation errors omitted.`);
    process.exitCode = 1;
} else console.log(`Workspace manifests, ${sourceOnly ? 'source aliases' : 'compiled exports'}, direct imports and package graph are valid (${packages.length} packages).`);
