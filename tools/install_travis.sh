#!/usr/bin/env bash

# Echo each command
set -x

# Exit on error.
set -e

if [[ "${AUDI_BUILD}" == manylinux* ]]; then
    cd ..;
    docker pull ${DOCKER_IMAGE};
    docker run --rm -e TWINE_PASSWORD -e AUDI_BUILD -e TRAVIS_TAG -v `pwd`:/audi $DOCKER_IMAGE bash /audi/tools/install_docker.sh
fi

set +e
set +x
