const common = require('./webpack.config.common.js');
const createApp = common.createApp;
module.exports = [
    createApp('playground', 'molstar')
];