/**
 * Copyright (c) 2017-2025 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Eric E <etongfu@@outlook.com>
 */
import * as esbuild from 'esbuild';
import * as fs from 'fs';
import * as path from 'path';
import * as argparse from 'argparse';
import { sassPlugin } from 'esbuild-sass-plugin';
import * as os from 'os';
import { execSync } from 'child_process';

const Apps = [
    // Apps
    { kind: 'app', name: 'viewer' },
    { kind: 'app', name: 'docking-viewer' },
    { kind: 'app', name: 'mesoscale-explorer' },
    { kind: 'app', name: 'mvs-stories', globalName: 'mvsStories', filename: 'mvs-stories.js' },

    // Examples
    { kind: 'example', name: 'proteopedia-wrapper' },
    { kind: 'example', name: 'basic-wrapper' },
    { kind: 'example', name: 'lighting' },
    { kind: 'example', name: 'alpha-orbitals' },
    { kind: 'example', name: 'alphafolddb-pae' },
    { kind: 'example', name: 'mvs-stories' },
    { kind: 'example', name: 'ihm-restraints' },
    { kind: 'example', name: 'interactions' },
    { kind: 'example', name: 'ligand-editor' },
];

function findApp(name, kind) {
    return Apps.find(a => a.name === name && a.kind === kind);
}

function mkDir(dir) {
    try {
        if (!fs.existsSync(dir)) {
            fs.mkdirSync(dir, { recursive: true });
        }
    } catch (error) {
        console.error(`Failed to create directory ${dir}:`, error);
        process.exit(1);
    }
}

function handleFileError(error, operation, path) {
    console.error(`Failed to ${operation} ${path}:`, error);
    process.exit(1);
}

function fileLoaderPlugin(options) {
    mkDir(options.out);

    return {
        name: 'file-loader',
        setup(build) {
            build.onLoad({ filter: /\.jpg$/ }, async (args) => {
                try {
                    const name = path.basename(args.path);
                    mkDir(path.resolve(options.out, 'images'));
                    await fs.promises.copyFile(args.path, path.resolve(options.out, 'images', name));
                    return {
                        contents: `images/${name}`,
                        loader: 'text',
                    };
                } catch (error) {
                    handleFileError(error, 'copy', args.path);
                }
            });
            build.onLoad({ filter: /\.(html|ico)$/ }, async (args) => {
                const name = path.basename(args.path);
                await fs.promises.copyFile(args.path, path.resolve(options.out, name));
                return {
                    contents: '',
                    loader: 'empty',
                };
            });
        },
    };
}

function examplesCssRenamePlugin({ root }) {
    return {
        name: 'example-css-rename',
        setup(build) {
            build.onEnd(async () => {
                if (fs.existsSync(path.resolve(root, 'index.css'))) {
                    await fs.promises.rename(
                        path.resolve(root, 'index.css'),
                        path.resolve(root, 'molstar.css')
                    );
                }
            });
        }
    };
}

async function createBundle(app) {
    const { name, kind } = app;

    const prefix = kind === 'app'
        ? `./build/${name}`
        : `./build/examples/${name}`;

    let entry = `./src/${kind}s/${name}/index.ts`;
    if (!fs.existsSync(entry)) {
        entry = `./src/${kind}s/${name}/index.tsx`;
    }

    let filename = app.filename;
    if (!filename) {
        filename = kind === 'app' ? 'molstar.js' : 'index.js';
    }

    const ctx = await esbuild.context({
        entryPoints: [entry],
        tsconfig: './tsconfig.json',
        bundle: true,
        minify: isProduction,
        minifyIdentifiers: false,
        sourcemap: includeSourceMap,
        globalName: app.globalName || 'molstar',
        outfile: kind === 'app'
            ? `./build/${name}/${filename}`
            : `./build/examples/${name}/${filename}`,
        plugins: [
            fileLoaderPlugin({ out: prefix }),
            sassPlugin({
                type: 'css',
                silenceDeprecations: ['import'],
                logger: {
                    warn: (msg) => console.warn(msg),
                    debug: () => { },
                }
            }),
            ...(kind === 'example' ? [examplesCssRenamePlugin({ root: prefix })] : []),
        ],
        external: ['crypto', 'fs', 'path', 'stream'],
        loader: {
        },
        color: true,
        logLevel: 'info',
        define: {
            'process.env.DEBUG': JSON.stringify(process.env.DEBUG || false),
            __MOLSTAR_PLUGIN_VERSION__: JSON.stringify(VERSION),
            __MOLSTAR_BUILD_TIMESTAMP__: `${TIMESTAMP}`,
        },
    });

    await ctx.rebuild();

    if (!isProduction) await ctx.watch();
}

const argParser = new argparse.ArgumentParser({
    add_help: true,
    description: 'Mol* Build'
});
argParser.add_argument('--prd', {
    help: 'Create a production build.',
    required: false,
    action: 'store_true',
});
argParser.add_argument('--src-map', {
    help: 'Include source map.',
    required: false,
    action: 'store_true',
});
argParser.add_argument('--apps', '-a', {
    help: 'Apps to build.',
    required: false,
    nargs: '*',
});
argParser.add_argument('--examples', '-e', {
    help: 'Examples to build.',
    required: false,
    nargs: '*',
});
argParser.add_argument('--port', '-p', {
    help: 'Port.',
    required: false,
    default: 1338,
    type: 'int',
});

argParser.add_argument('--host', {
    help: 'Show all available host addresses.',
    required: false,
    action: 'store_true',
});

const args = argParser.parse_args();


const isProduction = !!args.prd;
const includeSourceMap = !!args.src_map;

const VERSION = isProduction ? JSON.parse(fs.readFileSync('./package.json', 'utf8')).version : '(dev build)';
const TIMESTAMP = Date.now();

const apps = (!args.apps ? [] : (args.apps.length ? args.apps.map(a => findApp(a, 'app')).filter(a => a) : Apps.filter(a => a.kind === 'app')));
const examples = (!args.examples ? [] : (args.examples.length ? args.examples.map(e => findApp(e, 'example')).filter(a => a) : Apps.filter(a => a.kind === 'example')));

console.log('Apps:', apps.map(a => a.name));
console.log('Examples:', examples.map(e => e.name));
console.log('');

function getLocalIPs() {
    const interfaces = os.networkInterfaces();
    const ips = [];

    for (const name of Object.keys(interfaces)) {
        for (const iface of interfaces[name]) {
            // Skip internal and non-IPv4 addresses
            if (iface.internal || iface.family !== 'IPv4') continue;
            ips.push(iface.address);
        }
    }

    return ips;
}

function writeVersionFiles() {
    const file = `export var PLUGIN_VERSION = '${VERSION}';\nexport var PLUGIN_VERSION_DATE = new Date(${TIMESTAMP})`;
    const files = ['./lib/mol-plugin/version.js', './lib/commonjs/mol-plugin/version.js'];
    for (const f of files) {
        if (!fs.existsSync(f)) continue;
        fs.writeFileSync(f, file);
    }
}

async function main() {
    const promises = [];
    console.log(isProduction ? 'Building apps...' : 'Initial build...');

    for (const app of apps) promises.push(createBundle(app));
    for (const example of examples) promises.push(createBundle(example));

    await Promise.all(promises);

    if (isProduction) {
        console.log('Building library...');
        try {
            execSync('npm run build:lib', { stdio: 'inherit' });
            writeVersionFiles();
        } catch (error) {
            console.error('Failed to build library:');
            console.error(error);
            process.exit(1);
        }

        console.log('Done.');
        process.exit(0);
    }

    console.log('Initial build complete.');

    const ctx = await esbuild.context({});
    ctx.serve({
        servedir: './',
        port: args.port,
        host: '0.0.0.0', // Always listen on all interfaces
    });

    console.log('');
    console.log(`Server URL: http://localhost:${args.port}`);
    if (args.host) {
        console.log('Available host addresses:');
        const ips = getLocalIPs();
        ips.forEach(ip => console.log(`  http://${ip}:${args.port}`));
    }
    console.log('');
    console.log('Watching for changes...');
    console.log('');
    console.log('Press Ctrl+C to stop.');
}

main().catch((err) => {
    console.error(err);
    process.exit(1);
});