# Docker image for the multi-omics practicals

This build context creates an RStudio Server image for the Day 1 and Day 2 practicals of the course **Integrative Multi-Omic Data Analysis for Biological and Biomedical Data**.

## Build the image

From the folder containing this `Dockerfile`, run:

```bash
docker build -t multiomics-practical:2026 .
```

The build may take some time to download the image and install packages.

## Run RStudio Server

```bash
docker run --rm \
  -p 8787:8787 \
  -e PASSWORD=multiomics \
  multiomics-practical:2026
```

Then open:

```text
http://localhost:8787
```

Login credentials:

```text
Username: rstudio
Password: multiomics
```

## Optional: keep outputs after stopping the container

The datasets are already inside the image. To keep generated output files on your host machine, mount a local folder:

```bash
mkdir -p course_output

docker run --rm \
  -p 8787:8787 \
  -e PASSWORD=multiomics \
  -v "$(pwd)/course_output:/home/rstudio/course/output" \
  multiomics-practical:2026
```

Then save exported figures, reports, or model objects into `/home/rstudio/course/output`.

# Pre-built Docker images

For convenience, we provide pre-built Docker images on Docker Hub for both Apple Silicon Macs and AMD64/x86 systems.

## Apple Silicon Macs (M-series)

Pull the image:

```bash
docker pull lujunyan1118/multiomics-practical-mac:latest
```

Run the container:

```bash
docker run --rm -p 8787:8787 -e PASSWORD=multiomics \
lujunyan1118/multiomics-practical-mac:latest
```

---

## AMD64 / x86 Linux and Windows systems

Pull the image:

```bash
docker pull lujunyan1118/multiomics-practical-amd64:latest
```

Run the container:

```bash
docker run --rm -p 8787:8787 -e PASSWORD=multiomics \
lujunyan1118/multiomics-practical-amd64:latest
```

---
