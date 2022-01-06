/**
 * Copyright (c) 2021 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

export const PluginFeatureDetection = {
    get preferWebGl1() {
        if (typeof navigator === 'undefined' || typeof window === 'undefined') return false;

        // WebGL2 isn't working in MacOS 12.0.1 Safari 15.1, 15.2. It is working in Safari 15.4 tech preview, so disabling all versions before that.
        // prefer webgl 1 based on the userAgent substring
        const unpportedSafariVersions = [
            'Version/15.1 Safari',
            'Version/15.2 Safari',
            'Version/15.3 Safari'
        ];
        if (unpportedSafariVersions.some(v => navigator.userAgent.indexOf(v) > 0)) {
            return true;
        }

        // Check for iOS device which enabled WebGL2 recently but it doesn't seem
        // to be full up to speed yet.

        // adapted from https://stackoverflow.com/questions/9038625/detect-if-device-is-ios
        const isIOS = /iPad|iPhone|iPod/.test(navigator.userAgent);
        const isAppleDevice = navigator.userAgent.includes('Macintosh');
        const isTouchScreen = navigator.maxTouchPoints >= 4; // true for iOS 13 (and hopefully beyond)
        return !(window as any).MSStream && (isIOS || (isAppleDevice && isTouchScreen));
    },
    get wboit() {
        if (typeof navigator === 'undefined' || typeof window === 'undefined') return true;

        // disable Wboit in Safari 15
        return !/Version\/15.\d Safari/.test(navigator.userAgent);
    }
};