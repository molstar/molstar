import path from 'node:path';
import { fileURLToPath } from 'node:url';
import ts from 'typescript';

const root = path.resolve(path.dirname(fileURLToPath(import.meta.url)), '../..');
const config = path.resolve(process.argv[2] ?? path.join(root, 'tsconfig.json'));
const formatHost = {
    getCanonicalFileName: file => file,
    getCurrentDirectory: ts.sys.getCurrentDirectory,
    getNewLine: () => ts.sys.newLine,
};
const report = diagnostic => process.stderr.write(ts.formatDiagnosticsWithColorAndContext([diagnostic], formatHost));
const host = ts.createSolutionBuilderHost(ts.sys, undefined, report);
host.getParsedCommandLine = file => {
    const parsed = ts.getParsedCommandLineOfConfigFile(file, { skipLibCheck: false }, {
        ...ts.sys,
        onUnRecoverableConfigFileDiagnostic: report,
    });
    if (parsed?.options.tsBuildInfoFile) {
        parsed.options.tsBuildInfoFile = parsed.options.tsBuildInfoFile.replace(/\.tsbuildinfo$/, '.full.tsbuildinfo');
    }
    return parsed;
};
console.log('Checking all projects and imported declarations (skipLibCheck=false, forced rebuild).');
process.exitCode = ts.createSolutionBuilder(host, [config], { force: true }).build();
if (process.exitCode === 0) console.log('Full declaration check passed.');
