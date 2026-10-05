/** Compact and enumerate package export patterns using Node's subpath precedence. */
const escape = value => value.replace(/[.*+?^${}()|[\]\\]/g, '\\$&');
function capture(pattern, value) {
    const parts = pattern.split('*');
    if (parts.length === 1) return pattern === value ? '' : undefined;
    const match = value.match(new RegExp(`^${parts.map(escape).join('(.*)')}$`));
    return match && match.slice(1).every(part => part === match[1]) ? match[1] : undefined;
}
function substitute(value, matched) {
    if (typeof value === 'string') return value.replaceAll('*', matched);
    if (Array.isArray(value)) return value.map(item => substitute(item, matched));
    if (value && typeof value === 'object') return Object.fromEntries(Object.entries(value).map(([key, item]) => [key, substitute(item, matched)]));
    return value;
}
export function exportTargets(value, targets = []) {
    if (typeof value === 'string') targets.push(value);
    else if (value && typeof value === 'object') for (const item of Object.values(value)) exportTargets(item, targets);
    return targets;
}
export function resolveExport(exports, key) {
    if (Object.hasOwn(exports, key)) return exports[key];
    const patterns = Object.keys(exports).filter(pattern => pattern.includes('*')).sort((a, b) => b.indexOf('*') - a.indexOf('*') || b.length - a.length);
    for (const pattern of patterns) {
        const matched = capture(pattern, key);
        if (matched !== undefined) return substitute(exports[pattern], matched);
    }
}
export function expandExports(exports, files) {
    const keys = new Set(Object.keys(exports).filter(key => !key.includes('*')));
    for (const [pattern, value] of Object.entries(exports)) {
        if (!pattern.includes('*') || value === null) continue;
        if (pattern.split('*').length !== 2) throw new Error(`Export keys must contain a single wildcard: ${pattern}`);
        const matches = new Set();
        for (const target of exportTargets(value)) {
            if (!target.startsWith('./') || target.split('/').some(part => ['..', 'node_modules'].includes(part))) throw new Error(`Export target must stay package-local: ${target}`);
            for (const file of files) {
                const matched = capture(target, `./${file}`);
                if (matched !== undefined) matches.add(pattern.replace('*', matched));
            }
        }
        if (!matches.size) throw new Error(`Export pattern ${pattern} has no matching files`);
        for (const key of matches) keys.add(key);
    }
    return Object.fromEntries([...keys].sort().flatMap(key => {
        const value = resolveExport(exports, key);
        return value === null || value === undefined ? [] : [[key, value]];
    }));
}
function moduleExport(extension) {
    return { 'molstar-src': `./src/*.${extension}`, types: './lib/*.d.ts', import: './lib/*.js' };
}
export function compactExports(exports, sourceFiles) {
    if (Object.keys(exports).some(key => key.includes('*'))) return exports;
    const entries = Object.entries(exports);
    const canonical = (key, value, template, pattern = './*') => {
        const matched = capture(pattern, key);
        return matched !== undefined && JSON.stringify(value) === JSON.stringify(substitute(template, matched));
    };
    const counts = ['ts', 'tsx'].map(ext => [ext, entries.filter(([key, value]) => canonical(key, value, moduleExport(ext))).length]);
    const [extension, count] = counts.sort((a, b) => b[1] - a[1])[0];
    const result = {};
    const templates = [];
    if (count > 1) templates.push(['./*', moduleExport(extension)]);
    const sass = { 'molstar-src': './src/*.scss', sass: './src/*.scss', default: './lib/*.scss' };
    for (const [pattern, template] of [['./*.scss', sass], ['./*.css', './lib/*.css']]) {
        if (entries.filter(([key, value]) => canonical(key, value, template, pattern)).length > 1) templates.push([pattern, template]);
    }
    for (const [pattern, template] of templates) result[pattern] = template;
    for (const [key, value] of entries) {
        if (!templates.some(([pattern, template]) => canonical(key, value, template, pattern))) result[key] = value;
    }
    if (count > 1) {
        // Declaration sources must keep their explicit type-only entry, not a fake JS export.
        if (sourceFiles.some(file => file.endsWith('.d.ts'))) result['./*.d'] = null;
        for (const file of sourceFiles) {
            const parts = file.split('/'), test = parts.indexOf('_test');
            if (test >= 0) result[`./${parts.slice(0, test + 1).join('/')}/*`] = null;
            else if (/\.test\.tsx?$/.test(file)) result[`./${file.replace(/\.tsx?$/, '')}`] = null;
        }
    }
    return result;
}
