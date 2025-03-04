const webpack = require('webpack');
const path = require('path');
const packageJson = require('./package.json');

const serviceWorkerConfig = {
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
            }
        ]
    },
    plugins: [
        new webpack.DefinePlugin({
            VERSION: JSON.stringify(packageJson.version),
        }),
    ],
    resolve: {
        extensions: ['.ts', '.js']
    }
};

module.exports = [
    serviceWorkerConfig
];