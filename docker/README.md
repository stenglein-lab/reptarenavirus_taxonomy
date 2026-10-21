# Docker image needed to reproduce this paper

This directory contains a [dockerfile](./Dockerfile) that defines an image that can be 
used to reproduce this paper.  This image can be used by Docker or 
singularity via nextflow and will provide the tools needed to analyze
the data, generate figures, tables, and paper text, and render the paper.  

[This image](https://github.com/users/stenglein-lab/packages/container/package/aedes_isv_vt) provides the following tools:
- R and R packages, including:
  - tidyverse
  - knitr
  - patchwork
  - etc (see Dockerfile)
- Quarto 

This image contains well-defined versions of all the above tools. 

**It is not necessary to recreate this image; it already exists and is available via the Github Container Registry (GHCR), [here](https://github.com/users/stenglein-lab/packages/container/package/reptarenavirus_taxonomy).** 
This Dockerfile is included in this repository in order to document how the image was created.


## Steps to create and publish image

### Login to GHCR using a personal access token

Using a personal access token and actual github username

```
echo YOUR_PAT | docker login ghcr.io -u YOUR_GITHUB_USERNAME --password-stdin
```

### Build the image

```
docker buildx build --platform linux/amd64 -t ghcr.io/stenglein-lab/reptarenavirus_taxonomy:1.0 .
```

### Publish (push) the image to GHCR

```
docker push ghcr.io/stenglein-lab/reptarenavirus_taxonomy:1.0
```

After creating the package, it is necessary to change the visibility to public for others to be able to use it.  Do this via package settings on the [package's page](https://github.com/users/stenglein-lab/packages/container/reptarenavirus_taxonomy/settings) on github.
