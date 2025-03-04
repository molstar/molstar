import * as esbuild from 'esbuild';
import * as fs from 'fs';
import * as path from 'path';
import * as argparse from 'argparse';
import { sassPlugin } from 'esbuild-sass-plugin';

const AllApps = ['viewer', 'docking-viewer', 'mesoscale-explorer'];
const AllExamples = ['proteopedia-wrapper', 'basic-wrapper', 'lighting', 'alpha-orbitals', 'alphafolddb-pae', 'mvs-kinase-story', 'ihm-restraints'];

function mkDir(dir) {
    if (!fs.existsSync(dir)) {
        fs.mkdirSync(dir, { recursive: true });
    }
}

function fileLoaderPlugin(options) {
    mkDir(options.out, { recursive: true });

    return {
        name: 'file-loader',
        setup(build) {
            build.onLoad({ filter: /\.jpg$/ }, async (args) => {
                const name = path.basename(args.path);
                mkDir(path.resolve(options.out, 'images'));
                await fs.promises.copyFile(args.path, path.resolve(options.out, 'images', name));
                return {
                    contents: `images/${name}`,
                    loader: 'text',
                }
            });
            build.onLoad({ filter: /\.(html|ico)$/ }, async (args) => {
                const name = path.basename(args.path);
                await fs.promises.copyFile(args.path, path.resolve(options.out, name));
                return {
                    contents: '',
                    loader: 'empty',
                }
            });
        },
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

const args = argParser.parse_args();

const apps = (!args.apps ? [] : (args.apps.length ? args.apps : AllApps)).filter(a => AllApps.includes(a));
const examples = (!args.examples ? [] : (args.examples.length ? args.examples : AllExamples)).filter(e => AllExamples.includes(e));

console.log('Apps:', apps);
console.log('Examples:', examples);
console.log('');

const promises = [];
for (const app of apps) promises.push(watch(app, 'app'));
for (const example of examples) promises.push(watch(example, 'example'));

console.log('Initial build...');

await Promise.all(promises);
console.log('Done.');

const ctx = await esbuild.context({});
ctx.serve({
    servedir: './build',
    port: args.port,
});

console.log('');
console.log(`Serving on http://localhost:${args.port}`);
console.log('');
console.log('Watching for changes...');
console.log('');
console.log('Press Ctrl+C to stop.');