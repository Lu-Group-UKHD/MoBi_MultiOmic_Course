
# Docker Image for Day 3: Single-cell Integration

This directory contains the Docker environment for the Day 3: Single-cell Integration module.

Docker Hub Repository: `maxnuber/sc_integration_lecture`
Instructions


## 1. Prepare your local environment

Before running the container, ensure you have a directory to save your output files so they are not lost when the container stops.

`mkdir -p $(pwd)/course_output`


## 2. Run the Docker Container

To run the container with GPU acceleration (required for deep learning integration methods), use the following command:

```bash
docker run --rm -d \

--name sc_integration_course \

--gpus all \

-p 8787:8787 \

-p 8888:8888 \

-e PASSWORD=bioc \

-v "$(pwd)/course_output:/home/rstudio/project/output" \

maxnuber/sc_integration_lecture:latest
```


Command breakdown:

    --rm: Automatically remove the container when it exits.

    -d: Run in detached mode (in the background).

    --name: Assigns the name sc_integration_course to the container.

    --gpus all: Enables all available NVIDIA GPUs (requires NVIDIA Container Toolkit).

    -p 8787:8787: Maps the RStudio Server port.

    -p 8888:8888: Maps the Jupyter/Proxy port.

    -e PASSWORD=bioc: Sets the login password for RStudio to bioc.

    -v ...: Mounts your local folder to the container to persist your data.

## 3. Accessing the IDEs

    RStudio: Open your browser and go to http://localhost:8787

        Username: rstudio

        Password: bioc

    Jupyter: (If applicable) Open http://localhost:8888

## 4. Stopping the container

When you are finished, you can stop the container by running:

`docker stop sc_integration_course`