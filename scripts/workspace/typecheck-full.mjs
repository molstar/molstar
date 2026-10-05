import fs from 'node:fs';
import path from 'node:path';
import { randomUUID } from 'node:crypto';
import { spawnSync } from 'node:child_process';
import { fileURLToPath } from 'node:url';
import ts from '@typescript/typescript6';

const root = path.resolve(path.dirname(fileURLToPath(import.meta.url)), '../..');
const config = path.resolve(process.argv[2] ?? path.join(root, 'tsconfig.json'));
const compiler = path.join(root, 'node_modules/typescript/bin/tsc');
const suffix = randomUUID();
const projects = new Map();
const generated = [];
const formatHost = {
    getCanonicalFileName: file => file,
    getCurrentDirectory: ts.sys.getCurrentDirectory,
    getNewLine: () => ts.sys.newLine,
};
const report = diagnostic => process.stderr.write(ts.formatDiagnosticsWithColorAndContext([diagnostic], formatHost));

// Use the compatibility API only to read configs; TypeScript 7 performs every check and emit.
function collect(file) {
    if (fs.statSync(file).isDirectory()) file = path.join(file, 'tsconfig.json');
    file = path.resolve(file);
    if (projects.has(file)) return projects.get(file).full;
    const parsed = ts.getParsedCommandLineOfConfigFile(file, {}, {
        ...ts.sys,
        onUnRecoverableConfigFileDiagnostic: report,
    });
    if (!parsed || parsed.errors.length) {
        parsed?.errors.forEach(report);
        throw new Error(`Cannot read TypeScript project ${file}`);
    }
    // Keep configs beside the originals so default rootDir and @types lookup stay identical.
    const full = path.join(path.dirname(file), `.tsconfig.full-${suffix}.json`);
    const project = { full, parsed, references: [] };
    projects.set(file, project);
    project.references = (parsed.projectReferences ?? []).map(reference => ({ path: collect(reference.path) }));
    return full;
}

try {
    const fullConfig = collect(config);
    for (const [file, project] of projects) {
        const compilerOptions = { skipLibCheck: false };
        const buildInfo = ts.getTsBuildInfoEmitOutputFilePath(project.parsed.options);
        if (buildInfo) compilerOptions.tsBuildInfoFile = buildInfo.replace(/\.tsbuildinfo$/, '.full.tsbuildinfo');
        generated.push(project.full);
        fs.writeFileSync(project.full, JSON.stringify({
            extends: file,
            compilerOptions,
            references: project.references,
        }, null, 2) + '\n');
    }
    console.log('Checking all projects and imported declarations with TypeScript 7 (skipLibCheck=false, forced rebuild).');
    const result = spawnSync(process.execPath, [compiler, '-b', fullConfig, '--force'], { cwd: root, stdio: 'inherit' });
    if (result.error) throw result.error;
    process.exitCode = result.status ?? 1;
    if (process.exitCode === 0) console.log('Full declaration check passed.');
} finally {
    for (const file of generated) fs.rmSync(file, { force: true });
}
