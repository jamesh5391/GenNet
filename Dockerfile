# Use Miniconda as the base image
FROM continuumio/miniconda3

# Set the working directory
WORKDIR /workspace

# Copy your pip requirements file into the container
COPY requirements_GenNet.txt .

# Create the environment
RUN conda create -n env_GenNet python=3.10.12 -y

# Install pip packages
RUN /bin/bash -c "source activate env_GenNet && \
    pip install --upgrade pip && \
    pip install -r requirements_GenNet.txt"

# Copy your GenNet code into the image (optional)
COPY . /workspace

# Create input/output folders (in case not mounted)
RUN mkdir -p /workspace/processed_data /workspace/results

# Activate environment for future RUN/CMD
SHELL ["conda", "run", "-n", "env_GenNet", "/bin/bash", "-c"]

ENTRYPOINT ["conda", "run", "-n", "env_GenNet"]
CMD ["/bin/bash"]

