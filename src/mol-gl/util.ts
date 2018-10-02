/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

function debugTexture(imageData: ImageData, scale = 1) {
    const canvas = document.createElement('canvas')
    canvas.width = imageData.width
    canvas.height = imageData.height
    const ctx = canvas.getContext('2d')
    if (!ctx) throw new Error('Could not create canvas 2d context')
    ctx.putImageData(imageData, 0, 0)
    canvas.toBlob(imgBlob => {
        const objectURL = window.URL.createObjectURL(imgBlob)
        const img = document.createElement('img')
        img.src = objectURL
        img.style.width = imageData.width * scale + 'px'
        img.style.height = imageData.height * scale + 'px'
        img.style.position = 'absolute'
        img.style.top = '0px'
        img.style.left = '0px'
        document.body.appendChild(img)
    }, 'image/png')
}