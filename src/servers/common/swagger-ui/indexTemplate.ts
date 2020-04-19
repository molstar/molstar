export default `<!DOCTYPE html>
<html lang="en">
    <head>
        <meta charset="UTF-8">
        <title>\${title}</title>
        <link rel="stylesheet" type="text/css" href="\${apiPrefix}/swagger-ui.css" >
        \${shortcutIconLink}

        <style>
            html
            {
                box-sizing: border-box;
                overflow: -moz-scrollbars-vertical;
                overflow-y: scroll;
            }
            *,
            *:before,
            *:after
            {
                box-sizing: inherit;
            }
            body
            {
                margin:0;
                background: #fafafa;
            }
        </style>
    </head>

    <body>
        <div id="swagger-ui"></div>

        <script src="\${apiPrefix}/swagger-ui-bundle.js"> </script>
        <script src="\${apiPrefix}/swagger-ui-standalone-preset.js"> </script>
        <script>
            function HidePlugin() {
                // this plugin overrides some components to return nothing
                return {
                    components: {
                        Topbar: function () { return null },
                        Models: function () { return null },
                    }
                }
            }
            window.onload = function () {
                var ui = SwaggerUIBundle({
                    url: '\${openapiJsonUrl}',
                    validatorUrl: null,
                    docExpansion: 'list',
                    dom_id: '#swagger-ui',
                    deepLinking: true,
                    presets: [
                        SwaggerUIBundle.presets.apis,
                        SwaggerUIStandalonePreset
                    ],
                    plugins: [
                        SwaggerUIBundle.plugins.DownloadUrl,
                        HidePlugin
                    ],
                    layout: 'StandaloneLayout'
                })
                window.ui = ui
            }
        </script>
    </body>
</html>`;