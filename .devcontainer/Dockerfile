# See here for image contents: https://hub.docker.com/r/jupyter/datascience-notebook/

# Specify docker image from https://hub.docker.com/r/jupyter/datascience-notebook/tags.
# Avoid :latest tag. Tags that includes both linux/amd64 and linux/arm64/v8 are recommended.
FROM jupyter/datascience-notebook:2023-01-09

# We want to run common-debian.sh from here:
# https://github.com/microsoft/vscode-dev-containers/tree/main/script-library#development-container-scripts
# But that script assumes that the main non-root user (in this case jovyan)
# is in a group with the same name (in this case jovyan).  So we must first make that so.
COPY .devcontainer/library-scripts/common-debian.sh /tmp/library-scripts/
USER root
RUN apt-get update \
 && groupadd jovyan \
 && usermod -g jovyan -a -G users jovyan \
 && bash /tmp/library-scripts/common-debian.sh \
 && apt-get clean -y && rm -rf /var/lib/apt/lists/* /tmp/library-scripts

# Jupyter lab extensions
RUN mamba install -c conda-forge jupyterlab-git -y && \
    mamba install -c conda-forge jupytext==1.14.7 -y && \
    jupyter lab build

# Jupyter lab settings
# Alt + Space to send the highlighted code to console
# auto-close brackets + line numbers in notebooks and console
# turn off autosave
RUN printf ' { \
    "@jupyterlab/shortcuts-extension:shortcuts": {   \
        "shortcuts":[{"command": "notebook:run-in-console", "keys": ["Alt Space"], "selector": ".jp-Notebook.jp-mod-editMode"}, \
                     {"command": "runmenu:run", "keys": ["Alt Space"], "selector": "[data-jp-code-runner]"}] \
    }, \
    "@jupyterlab/notebook-extension:tracker": {  \
        "codeCellConfig": { "autoClosingBrackets": true, "lineNumbers": true }  \
    },   \
    "@jupyterlab/console-extension:tracker": {   \
        "promptCellConfig": { "autoClosingBrackets": true, "lineNumbers": true }   \
    },  \
    "@jupyterlab/docmanager-extension:plugin": {   \
        "autosave": false   \
    }  \
}' > /opt/conda/share/jupyter/lab/settings/overrides.json

# [Optional] If your Project.toml rarely change, uncomment this section to add them to the image.
ENV JULIA_PROJECT=.julia::/tmp/julia-tmp/
 COPY Project.toml /tmp/julia-tmp/
 COPY Manifest.toml /tmp/julia-tmp/
 RUN julia -e 'using Pkg; Pkg.activate("/tmp/julia-tmp"); Pkg.instantiate()'
#\
#    && rm -rf /tmp/julia-tmp

# R packages
RUN mamba install -c conda-forge r-ggpattern=1.0.1 r-igraph=1.4.1 \
    r-ggraph=2.1.0 r-ggmap=3.0.1 r-ggnewscale=0.4.8 r-maps=3.4.1 \
    r-geojsonio==0.11.3 r-sf==1.0-12 \
    -y
RUN mamba install -c conda-forge r-ggalluvial==0.12.5 -y
RUN mamba install -c r r-cmocean==0.3-1 -y
RUN pip install polars
# [Optional] If your pip requirements rarely change, uncomment this section to add them to the image.
# COPY requirements.txt /tmp/pip-tmp/
# RUN pip3 --disable-pip-version-check --no-cache-dir install -r /tmp/pip-tmp/requirements.txt \
#    && rm -rf /tmp/pip-tmp

# [Optional] Uncomment this section to install additional OS packages.
# RUN apt-get update && export DEBIAN_FRONTEND=noninteractive \
#     && apt-get -y install --no-install-recommends <your-package-list-here>
