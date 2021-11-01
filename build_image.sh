#! /bin/bash

# Remove docker and image if it was previously built
docker stop bsm_tutorial_docker
docker rm bsm_tutorial_docker
docker image rm bsm_tutorial

# Build the docker image
docker build --rm -t bsm_tutorial .

# Run the docker image as a container
# notebooks folder is an external volume - the notebook changes and 
# DataFiles produced from the notebooks are saved to your computer
docker run -it -p 8888:8888 --name='bsm_tutorial_docker' -e GRANT_SUDO=yes --user root -v $PWD/notebooks:/home/jovyan/bsm_tutorial/notebooks bsm_tutorial