# BSM Tutorial
These notebooks aim to provide a user-friendly way to understand borehole strainmeter data availability and processing. Jump to the last **Notebooks** section to see the details of each jupyter notebook. Below follows step-by-step instructions for downloading and installing necessary dependencies to run the notebooks in a Docker container. 

### 1. Install Docker
&ensp;&ensp; Pick your system and follow instructions on the Docker website. 
* **Mac** - https://docs.docker.com/docker-for-mac/install/ 
* **Windows** - https://docs.docker.com/docker-for-windows/install/ 
* **Ubuntu** - https://docs.docker.com/install/linux/docker-ce/ubuntu/ 

*Once installed, type `docker run hello-world` in terminal to check if installed correctly.*

More information on [getting started, testing your installation, and developing.](https://docs.docker.com/get-started/) 

For more information about **Docker**, check out our [short introduction](https://www.unavco.org/gitlab/GeoSciFramework/docker-tutorial/wikis/Introduction-to-Docker). 

### 2. Clone the Git Repository to your local machine 
`git clone https://gitlab.com/cehanagan/bsm_tutorial.git`

***

### **Build the Docker image**
In a terminal, navigate to the directory where the Dockerfile is stored (`bsm_tutorial`) and run:

`chmod +x build_image.sh`

`./build_image.sh`

This will use the Dockerfile to build the image with the jupyter notebook dependencies and spotl. When the image is run, the container is about ~5.7 GB.

This will also run the image and launch jupyter notebook. In the browser, type `localhost:8888` and enter the token `bsm`. Make sure no other notebooks are already running on this port prior to building the image. 

### **Run the container**
Once the image is built, if the container is not running, start it with the following command:

`docker run -it -p 8888:8888 --name='bsm_tutorial_docker' -e GRANT_SUDO=yes --user root -v ~/[your_host_directory]/notebooks:/home/jovyan/bsm_tutorial bsm_tutorial`

To access the container in another terminal, run:

`docker exec -it bsm_tutorial_docker /bin/bash`

### **Stop and remove the container and image**
To stop the container:

`docker stop bsm_tutorial_docker`

To remove the container:

`docker rm bsm_tutorial_docker`

To remove the image:

`docker image rm bsm_tutorial`

***

### **Notebooks**
The `notebooks/` directory contains four jupyter notebooks:
- NB1.ipynb: An introduction to borehole strainmeters, conventions, and available data - no code or output. 
- NB2.ipynb: Steps through picking a station to analyze and linearizing the data (produces "Level 1" data). 
- NB3.ipynb: Steps through calculating and applying barometric pressure, linear, and tidal corrections, and applying the regional orientation matrix (uses NB2 output files and produces 2 files - gauge strain and corrections and regional strain and corrections). 
- NB4.ipynb: Interactive plotting for NB3 output. 



