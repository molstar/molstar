const path = require('path');
const fs = require('fs');
const webpack = require('webpack');
const ExtraWatchWebpackPlugin = require('extra-watch-webpack-plugin');
const MiniCssExtractPlugin = require('mini-css-extract-plugin');
const VERSION = require('./package.json').version;

class VersionFilePlugin {
    apply() {
        fs.writeFileSync(
            path.resolve(__dirname, 'lib/mol-plugin/version.js'),
            `export var PLUGIN_VERSION = '${VERSION}';\nexport var PLUGIN_VERSION_DATE = new Date(typeof __MOLSTAR_DEBUG_TIMESTAMP__ !== 'undefined' ? __MOLSTAR_DEBUG_TIMESTAMP__ : ${new Date().valueOf()});`);
    }
}

const sharedConfig = {
    module: {
        rules: [
            {
                test: /\.(html|ico)$/,
                use: [{
                    loader: 'file-loader',
                    options: { name: '[name].[ext]' }
                }]
            },
            {
                test: /\.(s*)css$/,
                use: [
                    MiniCssExtractPlugin.loader,
                    { loader: 'css-loader', options: { sourceMap: false } },
                    { loader: 'sass-loader', options: { sourceMap: false } },
                ]
            },
            {
                test: /\.(jpg)$/i,
                type: 'asset/resource',
            },
        ]
    },
    plugins: [
        new ExtraWatchWebpackPlugin({
            files: [
                './lib/**/*.scss',
                './lib/**/*.html'
            ],
        }),
        new webpack.DefinePlugin({
            'process.env.DEBUG': JSON.stringify(process.env.DEBUG),
            '__MOLSTAR_DEBUG_TIMESTAMP__': webpack.DefinePlugin.runtimeValue(() => `${new Date().valueOf()}`, true)
        }),
        new MiniCssExtractPlugin({ filename: 'molstar.css' }),
        new VersionFilePlugin(),
    ],
    resolve: {
        modules: [
            'node_modules',
            path.resolve(__dirname, 'lib/')
        ],
        fallback: {
            fs: false,
            vm: false,
            buffer: false,
            crypto: require.resolve('crypto-browserify'),
            path: require.resolve('path-browserify'),
            stream: require.resolve('stream-browserify'),
        }
    },
    watchOptions: {
        aggregateTimeout: 750
    }
};

function createEntry(src, outFolder, outFilename, isNode) {
    return {
        target: isNode ? 'node' : void 0,
        entry: path.resolve(__dirname, `lib/${src}.js`),
        output: { filename: `${outFilename}.js`, path: path.resolve(__dirname, `build/${outFolder}`) },
        ...sharedConfig
    };
}

function createEntryPoint(name, dir, out, library) {
    return {
        entry: path.resolve(__dirname, `lib/${dir}/${name}.js`),
        output: { filename: `${library || name}.js`, path: path.resolve(__dirname, `build/${out}`), library: library || out, libraryTarget: 'umd', assetModuleFilename: 'images/[hash][ext][query]', 'publicPath': '' },
        ...sharedConfig
    };
}

function createNodeEntryPoint(name, dir, out) {
    return {
        target: 'node',
        entry: path.resolve(__dirname, `lib/${dir}/${name}.js`),
        output: { filename: `${name}.js`, path: path.resolve(__dirname, `build/${out}`) },
        externals: {
            argparse: 'require("argparse")',
            'node-fetch': 'require("node-fetch")',
            'util.promisify': 'require("util.promisify")',
            xhr2: 'require("xhr2")',
        },
        ...sharedConfig
    };
}

function createApp(name, library) { return createEntryPoint('index', `apps/${name}`, name, library); }
function createExample(name) { return createEntry(`examples/${name}/index`, `examples/${name}`, 'index'); }
function createBrowserTest(name) { return createEntryPoint(name, 'tests/browser', 'tests'); }
function createNodeApp(name) { return createNodeEntryPoint('index', `apps/${name}`, name); }

module.exports = {
    createApp,
    createEntry,
    createExample,
    createBrowserTest,
    createNodeEntryPoint,
    createNodeApp
};
