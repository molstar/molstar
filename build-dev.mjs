import * as esbuild from 'esbuild';
import * as fs from 'fs';
import * as path from 'path';
import * as argparse from 'argparse';
import { sassPlugin } from 'esbuild-sass-plugin';
import * as os from 'os';

const AllApps = [
    'viewer',
    'docking-viewer',
    'mesoscale-explorer'
];

const AllExamples = [
    'proteopedia-wrapper',
    'basic-wrapper',
    'lighting',
    'alpha-orbitals',
    'alphafolddb-pae',
    'mvs-kinase-story',
    'ihm-restraints',
    'interactions',
];

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
                try {
                    const name = path.basename(args.path);
                    await fs.promises.copyFile(args.path, path.resolve(options.out, name));
                    return {
                        contents: '',
                        loader: 'empty',
                    };
                } catch (error) {
                    handleFileError(error, 'copy', args.path);
                }
            });
        },
    };
}

function examplesCssRenamePlugin({ root }) {
    return {
        name: 'example-css-rename',
        setup(build) {
            build.onEnd(async () => {
                try {
                    const cssPath = path.resolve(root, 'index.css');
                    if (fs.existsSync(cssPath)) {
                        await fs.promises.rename(
                            cssPath,
                            path.resolve(root, 'molstar.css')
                        );
                    }
                } catch (error) {
                    handleFileError(error, 'rename', path.resolve(root, 'index.css'));
                }
            });
        }
    };
}

async function checkFileAccess(filePath, mode = fs.constants.R_OK) {
    try {
        await fs.promises.access(filePath, mode);
        return true;
    } catch (error) {
        console.error(`Cannot access file ${filePath}:`, error);
        return false;
    }
}

async function watch(name, kind) {
    const prefix = kind === 'app'
        ? `./build/${name}`
        : `./build/examples/${name}`;

    let entry = `./src/${kind}s/${name}/index.ts`;
    if (!fs.existsSync(entry)) {
        entry = `./src/${kind}s/${name}/index.tsx`;
    }

    if (!await checkFileAccess(entry)) {
        console.error(`Entry file not found or not accessible: ${entry}`);
        process.exit(1);
    }

    const outDir = path.dirname(kind === 'app'
        ? `./build/${name}/molstar.js`
        : `./build/examples/${name}/index.js`);

    try {
        await fs.promises.access(outDir, fs.constants.W_OK);
    } catch (error) {
        console.error(`Cannot write to output directory ${outDir}:`, error);
        process.exit(1);
    }

    const ctx = await esbuild.context({
        entryPoints: [entry],
        tsconfig: './tsconfig.json',
        bundle: true,
        globalName: 'molstar',
        outfile: kind === 'app'
            ? `./build/${name}/molstar.js`
            : `./build/examples/${name}/index.js`,
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
    return ctx;
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

const apps = (!args.apps ? [] : (args.apps.length ? args.apps : AllApps)).filter(a => AllApps.includes(a));
const examples = (!args.examples ? [] : (args.examples.length ? args.examples : AllExamples)).filter(e => AllExamples.includes(e));

console.log('Apps:', apps);
console.log('Examples:', examples);
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
    const contexts = [];
    let serverCtx;

    // Handle cleanup on process termination
    const cleanup = async () => {
        console.log('\nCleaning up...');
        for (const ctx of contexts) {
            try {
                await ctx.dispose();
            } catch (error) {
                console.error('Error during cleanup:', error);
            }
        }
        if (serverCtx) {
            try {
                await serverCtx.dispose();
            } catch (error) {
                console.error('Error during server cleanup:', error);
            }
        }
        process.exit(0);
    };

    process.on('SIGINT', cleanup);
    process.on('SIGTERM', cleanup);

    try {
        for (const app of apps) {
            const ctx = await watch(app, 'app');
            contexts.push(ctx);
        }
        for (const example of examples) {
            const ctx = await watch(example, 'example');
            contexts.push(ctx);
        }

        console.log('Initial build...');
        await Promise.all(promises);
        console.log('Done.');

        serverCtx = await esbuild.context({});
        serverCtx.serve({
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
    } catch (error) {
        console.error('Build failed:', error);
        await cleanup();
        process.exit(1);
    }
}

main().catch(error => {
    console.error('Build failed:', error);
    process.exit(1);
});