const { createApp, createExample, createBrowserTest } = require('./webpack.config.common.js');

const examples = ['proteopedia-wrapper', 'basic-wrapper', 'lighting', 'alpha-orbitals', 'alphafolddb-pae'];
const tests = [
    'font-atlas',
    'marching-cubes',
    'render-lines', 'render-mesh', 'render-shape', 'render-spheres', 'render-structure', 'render-structure-grid', 'render-text',
    'parse-xtc'
];

const path = require('path');

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
    },
    {
        entry: {
            sw: './src/apps/viewer/sw.ts'  // Separate entry for the service worker
        },
        output: {
            filename: '[name].js',  // No hash for the service worker
            path: path.resolve(__dirname, 'build/viewer')  // Output to 'build/viewer'
        },
        module: {
            rules: [
                {
                    test: /\.ts$/,
                    use: 'ts-loader',
                    exclude: '/node_modules'
                }
            ]
        },
        resolve: {
            extensions: ['.ts', '.js']
        },
        optimization: {
            minimize: false  // Do not minimize the service worker
        },
    }
];