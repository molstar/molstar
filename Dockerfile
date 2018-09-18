# This is to build a container that demos the Molstar Canvas prototype
# Source material: https://nodejs.org/en/docs/guides/nodejs-docker-webapp/
# Source material: https://derickbailey.com/2017/05/31/how-a-650mb-node-js-image-for-docker-uses-less-space-than-a-50mb-image/
# Source material: https://hub.docker.com/_/node/

# Use the slimed NodeJS source, yielding a space savings of 600MB (~66% of total)
FROM node:alpine

# Create app directory
WORKDIR /usr/src/app

# Install app dependencies
# A wildcard is used to ensure both package.json AND package-lock.json AND tslint.json AND tsconfig.json are copied
# where available (npm@5+)
COPY *.json ./

# Install all dependencies and copy results
RUN npm install
COPY . .

# Build application and bundle results
RUN npm run build
COPY build/ build/

# Build Canvas application and bundle results
RUN npm run build-canvas
COPY build/ build/

# Open ports for HTTP
EXPOSE 8080/tcp

# Setup standalone simple webserver to run the demo
RUN npm install http-server -g

# Start NodeJS at container stand up
CMD [ "http-server", "build/canvas/", "-p", "8080" ]

# Developer helpers (what is inside this container?)
RUN node -v
RUN ls -alh
