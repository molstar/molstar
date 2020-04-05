const createViewer = require('./webpack.config.viewer.js')[0];
module.exports = [
    {
        ...createViewer,
        devtool: 'eval'
    }
]