import fs from 'node:fs/promises';
import path from 'node:path';
import { spawnSync } from 'node:child_process';
import { fileURLToPath } from 'node:url';

const root = path.resolve(path.dirname(fileURLToPath(import.meta.url)), '../..');
const inventory = JSON.parse(await fs.readFile(path.join(root, 'scripts/workspace/inventory.json'), 'utf8'));
const version = JSON.parse(await fs.readFile(path.join(root, 'version.json'), 'utf8')).version;
const publicPackages = (inventory.packages ?? []).filter(pkg => pkg.private !== true);
const packageMap = new Map((inventory.packages ?? []).map(pkg => [pkg.name, pkg]));
const output = path.join(root, 'build/packages');
const requested = process.argv.slice(2).filter(arg => !arg.startsWith('--'));
const selected = requested.length ? publicPackages.filter(pkg => requested.includes(pkg.name) || requested.includes(pkg.path)) : publicPackages;
if (!selected.length) throw new Error(requested.length ? `No public packages match: ${requested.join(', ')}` : 'Inventory contains no public packages.');
await fs.rm(output, { recursive: true, force: true });
await fs.mkdir(output, { recursive: true });

function run(command, args, cwd) {
    const result = spawnSync(command, args, { cwd, encoding: 'utf8' });
    if (result.error) throw result.error;
    if (result.status !== 0) throw new Error(`${command} ${args.join(' ')} failed in ${cwd}:\n${result.stderr || result.stdout}`);
    return result.stdout;
}
function checkInternalRanges(manifest, pkg) {
    const errors = [];
    for (const field of ['dependencies', 'peerDependencies', 'optionalDependencies']) {
        for (const [name, range] of Object.entries(manifest[field] ?? {})) {
            if (packageMap.has(name) && range !== version) errors.push(`${pkg.name}: packed ${field}.${name} is ${range}, expected ${version}`);
        }
    }
    return errors;
}

function manifestTargets(manifest) {
    const targets = [];
    const visit = (value, conditions = [], acceptRelativeOnly = true) => {
        if (typeof value === 'string') {
            if (value.startsWith('./') || !acceptRelativeOnly) {
                targets.push({ path: value.replace(/^\.\//u, ''), conditions });
            }
            return;
        }
        if (!value || typeof value !== 'object') return;
        for (const [condition, target] of Object.entries(value)) visit(target, [...conditions, condition], acceptRelativeOnly);
    };
    for (const field of ['main', 'module', 'types', 'typings', 'style']) {
        const value = manifest[field];
        if (typeof value === 'string' && value.startsWith('./')) targets.push({ path: value.slice(2), conditions: [field] });
    }
    const bins = typeof manifest.bin === 'string' ? { [manifest.name]: manifest.bin } : manifest.bin;
    for (const [name, value] of Object.entries(bins ?? {})) {
        if (typeof value === 'string') targets.push({ path: value, conditions: [`bin:${name}`] });
    }
    for (const value of Object.values(manifest.exports ?? {})) visit(value, ['exports']);
    return targets;
}

function checkPackedContents(manifest, files, pkg) {
    const errors = [];
    const included = new Set(files);
    for (const target of manifestTargets(manifest)) {
        const normalized = path.posix.normalize(target.path.replace(/^\.\//u, ''));
        if (normalized.startsWith('../') || normalized === '..' || path.posix.isAbsolute(normalized)) {
            errors.push(`${pkg.name}: published target escapes package: ${target.path}`);
        } else if (!included.has(normalized)) {
            errors.push(`${pkg.name}: packed manifest target is missing from tarball: ${target.path} (${target.conditions.join('/')})`);
        }
    }
    for (const file of files) {
        const parts = file.split('/');
        const isTestTree = parts.includes('_test') || parts.includes('__tests__') || parts.includes('test') || parts.includes('tests') || parts.includes('fixtures');
        if (/\.tsbuildinfo$/iu.test(file)) errors.push(`${pkg.name}: tsbuildinfo unexpectedly packed: ${file}`);
        if (['src', 'lib'].includes(parts[0]) && (isTestTree || /(?:^|[/.])[^/]*\.test\.[^/]+$/iu.test(file))) {
            errors.push(`${pkg.name}: test fixture unexpectedly packed: ${file}`);
        }
    }
    return errors;
}

const failures = [];
for (const pkg of selected) {
    const packageDir = path.join(root, pkg.path);
    const manifestPath = path.join(packageDir, 'package.json');
    let sourceManifest;
    try { sourceManifest = JSON.parse(await fs.readFile(manifestPath, 'utf8')); } catch { failures.push(`${pkg.name}: missing or invalid ${pkg.path}/package.json`); continue; }
    if (sourceManifest.version !== version) failures.push(`${pkg.name}: source version ${sourceManifest.version} != ${version}`);
    const before = new Set((await fs.readdir(output)).filter(name => name.endsWith('.tgz')));
    run('pnpm', ['pack', '--pack-destination', output], packageDir);
    const archives = (await fs.readdir(output)).filter(name => name.endsWith('.tgz') && !before.has(name));
    if (archives.length !== 1) { failures.push(`${pkg.name}: expected one tarball, found ${archives.length}`); continue; }
    const archive = path.join(output, archives[0]);
    let packed;
    try { packed = JSON.parse(run('tar', ['-xOf', archive, 'package/package.json'], root)); } catch (error) { failures.push(`${pkg.name}: cannot inspect ${archives[0]}: ${error.message}`); continue; }
    let files;
    try {
        files = run('tar', ['-tzf', archive], root).split(/\r?\n/u).filter(Boolean).map(name => name.replace(/^package\//u, '').replace(/\/$/u, '')).filter(Boolean);
    } catch (error) { failures.push(`${pkg.name}: cannot list ${archives[0]}: ${error.message}`); continue; }
    if (packed.name !== pkg.name) failures.push(`${pkg.name}: packed manifest name is ${packed.name}`);
    if (packed.version !== version) failures.push(`${pkg.name}: packed version ${packed.version} != ${version}`);
    failures.push(...checkInternalRanges(packed, pkg));
    failures.push(...checkPackedContents(packed, files, pkg));
    console.log(`Packed ${pkg.name}@${packed.version} → ${path.relative(root, archive)}`);
}
if (failures.length) {
    console.error(failures.join('\n'));
    process.exitCode = 1;
} else console.log(`Verified ${selected.length} public package tarball(s) at version ${version}.`);
