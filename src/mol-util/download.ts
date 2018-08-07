/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

function openUrl (url: string) {
    const opened = window.open(url, '_blank')
    if (!opened) {
        window.location.href = url
    }
}

export function download (data: Blob|string, downloadName = 'download') {
    // using ideas from https://github.com/eligrey/FileSaver.js/blob/master/FileSaver.js

    if (!data) return

    const ua = window.navigator.userAgent
    const isSafari = /Safari/i.test(ua)
    const isChromeIos = /CriOS\/[\d]+/.test(ua)

    const a = document.createElement('a')

    function open (str: string) {
        openUrl(isChromeIos ? str : str.replace(/^data:[^;]*;/, 'data:attachment/file;'))
    }

    if (typeof navigator !== 'undefined' && navigator.msSaveOrOpenBlob) {
        // native saveAs in IE 10+
        navigator.msSaveOrOpenBlob(data, downloadName)
    } else if ((isSafari || isChromeIos) && FileReader) {
        if (data instanceof Blob) {
            // no downloading of blob urls in Safari
            const reader = new FileReader()
            reader.onloadend = () => open(reader.result as string)
            reader.readAsDataURL(data)
        } else {
            open(data)
        }
    } else {
        let objectUrlCreated = false
        if (data instanceof Blob) {
            data = URL.createObjectURL(data)
            objectUrlCreated = true
        }

        if ('download' in a) {
            // download link available
            a.style.display = 'hidden'
            document.body.appendChild(a)
            a.href = data
            a.download = downloadName
            a.target = '_blank'
            a.click()
            document.body.removeChild(a)
        } else {
            openUrl(data)
        }

        if (objectUrlCreated) {
            window.URL.revokeObjectURL(data)
        }
    }
}