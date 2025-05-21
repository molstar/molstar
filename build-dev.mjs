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

const Apps = [
    { name: 'viewer' },
    { name: 'docking-viewer' },
    { name: 'mesoscale-explorer' },
    { name: 'mvs-stories', globalName: 'mvsStories', filename: 'mvs-stories.js' },
];

const Examples = [
    { name: 'proteopedia-wrapper' },
    { name: 'basic-wrapper' },
    { name: 'lighting' },
    { name: 'alpha-orbitals' },
    { name: 'alphafolddb-pae' },
    { name: 'mvs-stories' },
    { name: 'ihm-restraints' },
    { name: 'interactions' },
    { name: 'ligand-editor' },
];

function findApp(name, kind) {
    if (kind === 'app') return Apps.find(a => a.name === name);
    if (kind === 'example') return Examples.find(e => e.name === name);
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

async function watch(app, kind) {
    const name = app.name;

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
    });

    await ctx.rebuild();
    await ctx.watch();
}

const argParser = new argparse.ArgumentParser({
    add_help: true,
    description: 'Mol* development build'
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

const apps = (!args.apps ? [] : (args.apps.length ? args.apps.filter(a => findApp(a, 'app')) : Apps));
const examples = (!args.examples ? [] : (args.examples.length ? args.examples.filter(e => findApp(e, 'example')) : Examples));

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

async function main() {
    const promises = [];
    for (const app of apps) promises.push(watch(app, 'app'));
    for (const example of examples) promises.push(watch(example, 'example'));

    console.log('Initial build...');

    await Promise.all(promises);
    console.log('Done.');

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

main().catch(console.error);