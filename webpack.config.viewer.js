const common = require('./webpack.config.common.js');
const createApp = common.createApp;
module.exports = [
    createApp('viewer', 'molstar'),
    createApp('mesoscale-explorer', 'molstar'),
];