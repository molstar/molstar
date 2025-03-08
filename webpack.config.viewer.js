const common = require('./webpack.config.common.js');
const VERSION = require('./package.json').version;
const createApp = common.createApp;
module.exports = (env, argv) => {
    return [
        createApp('viewer', 'molstar', [{
            from: 'lib/apps/viewer/sw.js',
            transform: (content) => {
                return content.toString().replace('const VERSION = \'\';', `const VERSION = '${VERSION}${env.WEBPACK_WATCH ? '-' + new Date().valueOf() : ''}';`);
            }
        }]),
        createApp('mesoscale-explorer', 'molstar'),
    ];
};