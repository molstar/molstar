const path = require('path');
const ExtraWatchWebpackPlugin = require('extra-watch-webpack-plugin');
const ExtractTextPlugin = require('extract-text-webpack-plugin');
module.exports = {
    module: {
        rules: [
            {
                loader: 'raw-loader',
                test: /\.(glsl|frag|vert)$/,
                include: [ path.resolve(__dirname, 'build/node_modules/') ],
            },
            {
                loader: 'glslify-loader',
                test: /\.(glsl|frag|vert)$/,
                include: [ path.resolve(__dirname, 'build/node_modules/') ]
            },

            {
                loader: 'file-loader',
                test: /\.(woff2?|ttf|otf|eot|svg|html)$/,
                include: [ path.resolve(__dirname, 'build/node_modules/') ],
                options: {
                    name: '[name].[ext]'
                }
            },
            {
                test:/\.(s*)css$/,
                use: ExtractTextPlugin.extract({
                    fallback:'style-loader',
                    use:['css-loader', 'resolve-url-loader', 'sass-loader'],
                })
            }
        ]
    },
    plugins: [
        new ExtraWatchWebpackPlugin({
            files: [
                './build/node_modules/**/*.vert',
                './build/node_modules/**/*.frag',
                './build/node_modules/**/*.glsl',
                './build/node_modules/**/*.scss',
                './build/node_modules/**/*.html'
            ],
        }),
        new ExtractTextPlugin({ filename:'app.css' }),
    ],
}