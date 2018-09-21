# This is to build a container that demos the Mol* canvas app
# Source material: https://nodejs.org/en/docs/guides/nodejs-docker-webapp/
# Source material: https://derickbailey.com/2017/05/31/how-a-650mb-node-js-image-for-docker-uses-less-space-than-a-50mb-image/
# Source material: https://hub.docker.com/_/node/


### Builder Stage (Multi-Stage Docker Build)
# Use the slimed NodeJS source, yielding a space savings of 600MB (~66% of total)
FROM node:alpine as builder
ENV NODEROOT /usr/src/app/

# Create app directory
WORKDIR $NODEROOT

# Install app dependencies
# A wildcard is used to ensure the following are copied
# package.json AND package-lock.json AND tslint.json AND tsconfig.json are copied where available (npm@5+)
COPY *.json ./

# Install all dependencies and copy results
RUN npm install
COPY . .

# Build library and canvas application then copy results
RUN npm run build
RUN npm run build-canvas
COPY . .

### Runtime Stage (Multi-Stage Docker Build)
FROM httpd:2.4 as runtime
ENV NODEROOT /usr/src/app
ENV HTTPDROOT /usr/local/apache2/htdocs/

## Code must be placed into /usr/local/apache2/htdocs/
WORKDIR $HTTPDROOT
COPY --from=builder $NODEROOT/build/canvas/ .

# Open ports for HTTP
EXPOSE 80
