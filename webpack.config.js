const { createApp, createExample, createBrowserTest } = require('./webpack.config.common.js');

const examples = ['proteopedia-wrapper', 'basic-wrapper', 'lighting', 'alpha-orbitals', 'alphafolddb-pae'];
const tests = [
    'font-atlas',
    'marching-cubes',
    'render-lines', 'render-mesh', 'render-shape', 'render-spheres', 'render-structure', 'render-structure-grid', 'render-text',
    'parse-xtc'
];

const path = require('path');
const webpack = require('webpack');
const packageJson = require('./package.json');

function createApp(name, library) {
    return {
        entry: {
            main: `./src/apps/${name}/index.ts`,
            sw: `./src/apps/${name}/sw.ts`  // Ensure service worker is included
        },
        output: {
            filename: '[name].[contenthash].js',  // Include hash for cache invalidation
            path: path.resolve(__dirname, `build/${name}`)
        },
        module: {
            rules: [
                {
                    test: /\.ts$/,
                    use: 'ts-loader',
                    exclude: /node_modules/
                }
            ]
        },
        resolve: {
            extensions: ['.ts', '.js']
        },
        plugins: [
            new webpack.DefinePlugin({
                'process.env.VERSION': JSON.stringify(packageJson.version),
            }),
        ],
        optimization: {
            minimize: true
        }
    };
}

module.exports = [
    createApp('viewer', 'molstar'),
    createApp('docking-viewer', 'molstar'),
    createApp('mesoscale-explorer', 'molstar'),
    ...examples.map(createExample),
    ...tests.map(createBrowserTest),
    {
        entry: {
            main: './src/index.ts'
        },
        output: {
            filename: '[name].[contenthash].js',  // Include hash for cache invalidation
            path: path.resolve(__dirname, 'build/viewer')
        },
        module: {
            rules: [
                {
                    test: /\.ts$/,
                    use: 'ts-loader',
                    exclude: /node_modules/
                }
            ]
        },
        resolve: {
            extensions: ['.ts', '.js']
        },
        optimization: {
            minimize: true
        },
        plugins: [
            new webpack.DefinePlugin({
                'process.env.VERSION': JSON.stringify(packageJson.version),
            }),
        ]
    },
    {
        entry: {
            sw: './src/apps/viewer/sw.ts'  // Entry point for the service worker
        },
        output: {
            filename: '[name].js',
            path: path.resolve(__dirname, 'build/viewer')
        },
        module: {
            rules: [
                {
                    test: /\.ts$/,
                    use: 'ts-loader',
                    exclude: /node_modules/
                },
                {
                    test: /\.json$/,
                    type: 'json'
                }
            ]
        },
        resolve: {
            extensions: ['.ts', '.js', '.json']
        },
        optimization: {
            minimize: false  // Do not minimize the service worker
        },
        plugins: [
            new webpack.DefinePlugin({
                'process.env.VERSION': JSON.stringify(packageJson.version),
            }),
        ]
    }
];